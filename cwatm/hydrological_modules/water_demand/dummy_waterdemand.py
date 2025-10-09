    def dynamic(self):
        """
        Dynamic part of the water demand module

        * calculate the fraction of water from surface water vs. groundwater
        * get non-Irrigation water demand and its return flow fraction
        """

        if checkOption('includeWaterDemand'):

            # for debugging of a specific date
            # if (globals.dateVar['curr'] >= 137):
            #    ii =1

            # ----------------------------------------------------
            # WATER DEMAND

            # Fix year of water demand on predefined year
            wd_date = globals.dateVar['currDate']
            if self.var.waterdemandFixed:
                wd_date = wd_date.replace(day=1)
                wd_date = wd_date.replace(year=self.var.waterdemandFixedYear)

            # calculate water demand for different sectors
            self.irrigation.dynamic()
            self.environmental_need.dynamic()
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.domestic.dynamic(wd_date)
                self.industry.dynamic(wd_date)
                self.livestock.dynamic(wd_date)

            # calculate total water demand
            talDemand, frac_industry, frac_domestic, frac_livestock = self.calc_totalDemand()

            # ----------------------------------------------------
            # WATER AVAILABILITY

            # conversion m3 -> m # minus environmental flow
            self.var.readAvlChannelStorageM = np.maximum(0.,
                                                         self.var.channelStorage * self.var.M3toM - self.var.envFlow)  # in [m]

            # -------------------------------------
            # WATER DEMAND vs. WATER AVAILABILITY
            # -------------------------------------

            if dateVar['newStart'] or dateVar['newMonth']:
                self.var.act_irrWithdrawalSW_month = globals.inZero.copy()
                self.var.act_irrWithdrawalGW_month = globals.inZero.copy()

                if 'commandAreasRelaxGwAbstraction' in binding and self.var.sectorSourceAbstractionFractions:
                    self.var.gwAbstractionFraction_Irrigation = np.where(self.var.reservoir_command_areas > 0,
                                                                         0.01,
                                                                         self.var.gwAbstractionFraction_Irrigation)

            if self.var.sectorSourceAbstractionFractions and 'commandAreasRelaxGwAbstraction' in binding and \
                    self.var.using_reservoir_command_areas:

                if dateVar['currDate'].day > 15:
                    self.var.gwAbstractionFraction_Irrigation = np.where(self.var.reservoir_command_areas > 0,
                                                                         loadmap('commandAreasRelaxGwAbstraction'),
                                                                         self.var.gwAbstractionFraction_Irrigation)
            
            
            # Desalination
            self.calc_desalination()
                
            
            # surface water abstraction that can be extracted to fulfill totalDemand
            # - based on ChannelStorage and swAbstractionFraction * totalDemand
            # sum up potential surface water abstraction (no groundwater abstraction under water and sealed area)

            # calculate channel abstraction for given sector-specific source abstraction fractions
            self.calc_channelAbstract_SourceAbstraction()

            # channel abstraction from lift areas
            self.calc_channelAbstract_Lift()


            if checkOption('includeWaterBodies'):
                self.var.abstractedLakeReservoirM3C = np.compress(self.var.compress_LR, globals.inZero.copy())                        
                        
                # water that is still needed from surface water
                remainNeed0 = np.maximum(pot_SurfaceAbstract - self.var.act_SurfaceWaterAbstract, 0)
                
                self.var.abstractedLakeReservoirM3C = np.compress(self.var.compress_LR, globals.inZero.copy())  
                mskWtrBody_unrestricted = self.var.waterBodyBuffer > 0
                
                # first from big Lakes and reservoirs, big lakes cover several gridcells
                # collect all water demand from lake pixels of the same id

                # remainNeedBig = npareatotal(remainNeed, self.var.waterBodyID)
                # not only the lakes and reservoirs but the command areas around water bodies e.g. here a buffer
                remainNeedBig = npareatotal(remainNeed0, self.var.waterBodyBuffer)
                remainNeedBigC = np.compress(self.var.compress_LR, remainNeedBig)

                # Storage of a big lake
                lakeResStorageC = np.where(self.var.waterBodyTypCTemp == 0, 0.,
                                           np.where(self.var.waterBodyTypCTemp == 1, self.var.lakeStorageC,
                                                    self.var.reservoirStorageM3C)) / self.var.MtoM3C

                minlake = np.maximum(0., 0.98 * lakeResStorageC)  # reasonable but arbitrary limit
                act_bigLakeAbstC = np.minimum(minlake, remainNeedBigC)
                # substract from both, because it is sorted by self.var.waterBodyTypCTemp
                self.var.lakeStorageC = self.var.lakeStorageC - act_bigLakeAbstC * self.var.MtoM3C
                self.var.lakeVolumeM3C = self.var.lakeVolumeM3C - act_bigLakeAbstC * self.var.MtoM3C
                self.var.reservoirStorageM3C = self.var.reservoirStorageM3C - act_bigLakeAbstC * self.var.MtoM3C
                # and from the combined one for waterbalance issues
                self.var.lakeResStorageC = self.var.lakeResStorageC - act_bigLakeAbstC * self.var.MtoM3C

                self.var.abstractedLakeReservoirM3C = act_bigLakeAbstC.copy() * self.var.MtoM3C

                self.var.lakeResStorage = globals.inZero.copy()
                np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)
                bigLakesFactorC = divideValues(act_bigLakeAbstC, remainNeedBigC)

                # and back to the big array
                bigLakesFactor = globals.inZero.copy()
                np.put(bigLakesFactor, self.var.decompress_LR, bigLakesFactorC)

                # bigLakesFactorAllaroundlake = npareamaximum(bigLakesFactor, self.var.waterBodyID)
                bigLakesFactorAllaroundlake = npareamaximum(bigLakesFactor, self.var.waterBodyBuffer)

                # abstraction from big lakes is partioned to the users around the lake
                self.var.act_bigLakeResAbst = remainNeed0  * mskWtrBody_unrestricted * bigLakesFactorAllaroundlake   

                # remaining need is used from small lakes
                remainNeed1 = remainNeed0 * (1 - bigLakesFactorAllaroundlake)
                # minlake = np.maximum(0.,self.var.smalllakeStorage - self.var.minsmalllakeStorage) * self.var.M3toM

                if returnBool('useSmallLakes'):
                    minlake = np.maximum(0., 0.98 * self.var.smalllakeStorage) * self.var.M3toM
                    self.var.act_smallLakeResAbst = np.minimum(minlake, remainNeed1)
                    # act_smallLakesres is substracted from small lakes storage
                    self.var.smalllakeVolumeM3 = self.var.smalllakeVolumeM3 - self.var.act_smallLakeResAbst * self.var.MtoM3
                    self.var.smalllakeStorage = self.var.smalllakeStorage - self.var.act_smallLakeResAbst * self.var.MtoM3
                else:
                    self.var.act_smallLakeResAbst = 0

                # available surface water is from river network + large/small lake & reservoirs
                self.var.act_SurfaceWaterAbstract = self.var.act_SurfaceWaterAbstract + self.var.act_bigLakeResAbst \
                                                    + self.var.act_smallLakeResAbst
                self.var.act_lakeAbst = self.var.act_bigLakeResAbst + self.var.act_smallLakeResAbst

                # Transfer water between reservoirs
                # Send storage between reservoirs using the Excel sheet reservoir_transfers within cwatm_settings.xlsx
                # Using the waterBodyIDs defined in the settings, designate
                # the Giver, the Receiver, and the daily fraction of live storage the Giver sends to the Receiver.
                # If the Receiver is already at capacity, the Giver does not send any storage.
                # Reservoirs can only send to one reservoir. Reservoirs can receive from several reservoirs.

                if 'reservoir_transfers' in option:
                    if checkOption('reservoir_transfers'):

                        for transfer in self.var.reservoir_transfers:

                            self.var.inZero_C = np.compress(self.var.compress_LR, globals.inZero.copy())

                            if returnBool('dynamicLakesRes'):
                                year = dateVar['currDate'].year
                            else:
                                year = loadmap('fixLakesResYear')

                            if transfer[0] > 0:
                                index_giver = np.where(self.var.waterBodyID_C == transfer[0])[0][0]
                                giver_already_constructed = self.var.resYearC[index_giver] <= year
                            else:
                                giver_already_constructed = True

                            if transfer[1] > 0:
                                index_receiver = np.where(self.var.waterBodyID_C == transfer[1])[0][0]
                                receiver_already_constructed = self.var.resYearC[index_receiver] <= year
                            else:
                                receiver_already_constructed = True

                            if receiver_already_constructed and giver_already_constructed:

                                reservoir_unused = self.var.resVolumeC - self.var.reservoirStorageM3C
                                if transfer[1] > 0:
                                    reservoir_unused_receiver = reservoir_unused[index_receiver]
                                else:
                                    reservoir_unused_receiver = 10e12

                                if transfer[0] == 0:
                                    # In this case, the fraction refers to the fraction of the receiver,
                                    # as the giver is infinite
                                    reservoir_storage_giver = self.var.resVolumeC[index_receiver]
                                else:
                                    reservoir_storage_giver = self.var.reservoirStorageM3C[index_giver]

                                reservoir_transfer_actual = np.minimum(reservoir_unused_receiver * 0.95,
                                                                       np.where(transfer[2] <= 1,
                                                                                reservoir_storage_giver * transfer[2],
                                                                                transfer[2]))

                                # print(transfer[0], 'donated', reservoir_transfer_actual, 'm3 to', transfer[1])

                                if transfer[0] > 0:  # There is a giver, not the ocean
                                    self.var.inZero_C[index_giver] = -reservoir_transfer_actual  # giver
                                    self.var.reservoir_transfers_out_M3C[index_giver] += reservoir_transfer_actual
                                else:
                                    self.var.reservoir_transfers_from_outside_M3C[index_receiver] \
                                        += reservoir_transfer_actual

                                if transfer[1] > 0:  # There is a receiver, not the ocean
                                    self.var.inZero_C[index_receiver] = reservoir_transfer_actual  # receiver
                                    self.var.reservoir_transfers_in_M3C[index_receiver] += reservoir_transfer_actual
                                else:
                                    self.var.reservoir_transfers_to_outside_M3C[index_giver] \
                                        += reservoir_transfer_actual



                                self.var.lakeStorageC += self.var.inZero_C
                                self.var.lakeVolumeM3C += self.var.inZero_C
                                self.var.lakeResStorageC += self.var.inZero_C
                                self.var.reservoirStorageM3C += self.var.inZero_C

                                self.var.reservoir_transfers_net_M3C += self.var.inZero_C
                                # Cancels out positive and negative if both receiving and giving

                                if transfer[1] == 0:
                                    to_outside_basin = globals.inZero.copy()
                                    np.put(to_outside_basin, self.var.decompress_LR, self.var.inZero_C)
                                    pot_Lake_Industry -= to_outside_basin * self.var.M3toM
                                    # self.var.Lake_Industry is updated below
                                    self.var.act_lakeAbst -= to_outside_basin * self.var.M3toM

                                    self.var.act_SurfaceWaterAbstract -= to_outside_basin * self.var.M3toM
                                    self.var.act_bigLakeResAbst -= to_outside_basin * self.var.M3toM
                                    
                                    self.var.industryDemand -= to_outside_basin * self.var.M3toM
                                    self.var.pot_industryConsumption -= to_outside_basin * self.var.M3toM
                                    self.var.ind_efficiency = divideValues(self.var.pot_industryConsumption,
                                                                           self.var.industryDemand)


                        np.put(self.var.reservoir_transfers_net_M3, self.var.decompress_LR,
                               self.var.reservoir_transfers_net_M3C)
                        self.var.reservoir_transfers_net_M3C = np.compress(self.var.compress_LR,
                                                                           globals.inZero.copy())

                        np.put(self.var.reservoir_transfers_in_M3, self.var.decompress_LR,
                               self.var.reservoir_transfers_in_M3C)
                        self.var.reservoir_transfers_in_M3C = np.compress(self.var.compress_LR,
                                                                          globals.inZero.copy())

                        np.put(self.var.reservoir_transfers_out_M3, self.var.decompress_LR,
                               self.var.reservoir_transfers_out_M3C)
                        self.var.reservoir_transfers_out_M3C = np.compress(self.var.compress_LR,
                                                                           globals.inZero.copy())

                        np.put(self.var.reservoir_transfers_from_outside_M3, self.var.decompress_LR,
                               self.var.reservoir_transfers_from_outside_M3C)
                        self.var.reservoir_transfers_from_outside_M3C = np.compress(self.var.compress_LR,
                                                                                    globals.inZero.copy())

                        np.put(self.var.reservoir_transfers_to_outside_M3, self.var.decompress_LR,
                               self.var.reservoir_transfers_to_outside_M3C)
                        self.var.reservoir_transfers_to_outside_M3C = np.compress(self.var.compress_LR,
                                                                                  globals.inZero.copy())

                        ###

                        if self.var.sectorSourceAbstractionFractions:

                            self.var.swAbstractionFraction_Res_Industry = \
                                np.where(self.var.reservoir_transfers_to_outside_M3 > 0, 0,
                                         self.var.swAbstractionFraction_Res_Industry)
                            self.var.gwAbstractionFraction_Industry = \
                                np.where(self.var.reservoir_transfers_to_outside_M3 > 0, 0,
                                         self.var.gwAbstractionFraction_Industry)
                        else:
                            pot_SurfaceAbstract -= to_outside_basin
                            # to avoid groundwater abstraction
                            self.var.swAbstractionFraction = \
                                np.where(self.var.reservoir_transfers_to_outside_M3 != 0, 1,
                                         self.var.swAbstractionFraction_nonIrr)

                # -------------------------------------

                if self.var.sectorSourceAbstractionFractions:

                    # A
                    self.var.Lake_Domestic = np.minimum(self.var.act_lakeAbst, pot_Lake_Domestic)
                    self.var.Lake_Livestock = np.minimum(self.var.act_lakeAbst - self.var.Lake_Domestic,
                                                         pot_Lake_Livestock)
                    self.var.Lake_Industry = np.minimum(
                        self.var.act_lakeAbst - self.var.Lake_Domestic - self.var.Lake_Livestock,
                        pot_Lake_Industry)
                    self.var.Lake_Irrigation = np.minimum(
                        self.var.act_lakeAbst - self.var.Lake_Domestic - self.var.Lake_Livestock - self.var.Lake_Industry,
                        pot_Lake_Irrigation)

                    # B
                    pot_Res_Domestic = np.minimum(
                        self.var.swAbstractionFraction_Res_Domestic * self.var.domesticDemand,
                        self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic - \
                            self.var.Lift_Domestic - self.var.Lake_Domestic)

                    pot_Res_Livestock = np.minimum(
                        self.var.swAbstractionFraction_Res_Livestock * self.var.livestockDemand,
                        self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock - \
                            self.var.Lift_Livestock - self.var.Lake_Livestock)

                    pot_Res_Industry = np.minimum(
                        self.var.swAbstractionFraction_Res_Industry * self.var.industryDemand,
                        self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry - \
                            self.var.Lift_Industry - self.var.Lake_Industry)

                    pot_Res_Irrigation = np.minimum(
                        self.var.swAbstractionFraction_Res_Irrigation * self.var.totalIrrDemand,
                        self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation - \
                            self.var.Lift_Irrigation - self.var.Lake_Irrigation)

                    # remainNeed2 = pot_Res_Domestic + pot_Res_Livestock + pot_Res_Industry + pot_Res_Irrigation
                    remainNeed2 = pot_Res_Irrigation

                else:

                    remainNeed2 = pot_SurfaceAbstract - self.var.act_SurfaceWaterAbstract

                #if self.var.load_command_areas:

                # The command area of a reservoir is the area that can receive water from this reservoir. Before
                # this step, each cell has attempted to satisfy its demands with local water using in-cell
                # channel, lift area, and lake. The remaining demand within each command area is totaled and
                # requested from the associated reservoir. The reservoir offers this water up to a daily maximum
                # relating to the available storage in the reservoir, defined in the Reservoir_releases_input_file.
                #
                # SETTINGS FILE AND INPUTS

                # -Activating
                # In the OPTIONS section towards the beginning of the settings file, add/set
                # using_reservoir_command_areas = True

                # - Command areas raster map Anywhere after the OPTIONS section (in WATERDEMAND, for example),
                # add/set reservoir_command_areas to a path holding... information about the command areas. This
                # Command areas raster map should assign the same positive integer coding to each cell within the
                # same segment. All other cells must Nan values, or values <= 0.

                # -Optional inputs
                #
                # Anywhere after the OPTIONS section, add/set Reservoir_releases_input_file to a path holding
                # information about irrigation releases. This should be a raster map (netCDF) of 366 values
                # determining the maximum fraction of available storage to be used for meeting water demand... in
                # the associated command area on the day of the year. If this is not included, a value of 0.01
                # will be assumed, i.e. 1% of the reservoir storage can be at most released into the command area
                # on each day.
                self.var.act_ResAbst = globals.inZero.copy()
                if self.var.sectorSourceAbstractionFractions:

                    # Domestic, livestock, and industrial demands are satisfied before irrigation

                    remainNeedPre = pot_Res_Domestic + pot_Res_Livestock + pot_Res_Industry
                    #print('water_demand.py: np.sum(remainNeedPre) with reservoirs', np.sum(remainNeedPre))

                    demand_Segment = np.where(self.var.reservoir_command_areas > 0,
                                              npareatotal(remainNeedPre * self.var.cellArea,
                                                          self.var.reservoir_command_areas),
                                              0)  # [M3]

                    #print('water_demand.py: np.sum(demand_Segment) with reservoirs', np.sum(demand_Segment))

                    # Reservoir associated with the Command Area
                    #
                    # If there is more than one reservoir in a command area, the storage of the reservoir with
                    # maximum storage in this time-step is chosen. The map resStorageTotal_alloc holds this
                    # maximum reservoir storage within a command area in all cells within that command area

                    # filter reservoirs so only non-restricted res. are acccounted for
                    ReservoirsThatAreCurrentlyReservoirs = np.where(self.var.waterBodyTypCTemp == 2, \
                                self.var.reservoirStorageM3C, np.where(self.var.waterBodyTypCTemp == 4, self.var.reservoirStorageM3C, 0))
                    ReservoirsThatAreCurrentlyReservoirs = np.where(np.compress(self.var.compress_LR, self.var.resId_restricted) == 0, \
                        ReservoirsThatAreCurrentlyReservoirs, 0)

                    reservoirStorageM3 = globals.inZero.copy()
                    # np.put(reservoirStorageM3, self.var.decompress_LR, self.var.reservoirStorageM3C)
                    np.put(reservoirStorageM3, self.var.decompress_LR, ReservoirsThatAreCurrentlyReservoirs)

                    resStorageTotal_alloc = np.where(self.var.reservoir_command_areas > 0,
                                                     npareamaximum(reservoirStorageM3,
                                                                   self.var.reservoir_command_areas), 0)  # [M3]

                    # In the map resStorageTotal_allocC, the maximum storage from each allocation segment is held
                    # in all reservoir cells within that allocation segment. We now correct to remove the
                    # reservoirs that are not this maximum-storage-reservoir for the command area.
                    resStorageTotal_allocC = np.compress(self.var.compress_LR, resStorageTotal_alloc)
                    resStorageTotal_allocC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C,
                                                         resStorageTotal_allocC)

                    day_of_year = globals.dateVar['currDate'].timetuple().tm_yday


                    if 'Reservoir_releases' in binding:
                        # resStorage_maxFracForIrrigation = 0.5 + globals.inZero.copy()
                        resStorage_maxFracForIrrigation = readnetcdf2('Reservoir_releases', day_of_year,
                                                                      useDaily='DOY', value='Fraction of Volume')
                        resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR,
                                                                       resStorage_maxFracForIrrigation)
                    elif self.var.reservoir_releases_excel_option:
                        resStorage_maxFracForIrrigation = globals.inZero.copy()
                        resStorage_maxFracForIrrigationC = np.where(self.var.lakeResStorage_release_ratioC > -1,
                                                                    self.var.reservoir_supply[dateVar['doy']-1],
                                                                    0.03)
                    else:
                        resStorage_maxFracForIrrigation = 0.03 + globals.inZero.copy()
                        resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR,
                                                                       resStorage_maxFracForIrrigation)

                    # resStorage_maxFracForIrrigationC holds the fractional rules found for each reservoir,
                    # so we must null those that are not the maximum-storage reservoirs

                    resStorage_maxFracForIrrigationC = np.multiply(
                        resStorageTotal_allocC == self.var.reservoirStorageM3C, resStorage_maxFracForIrrigationC)

                    np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR,
                           resStorage_maxFracForIrrigationC)


                    resStorage_maxFracForIrrigation_CA = np.where(self.var.reservoir_command_areas > 0,
                                                                  npareamaximum(resStorage_maxFracForIrrigation,
                                                                                self.var.reservoir_command_areas),
                                                                  0)

                    act_bigLakeResAbst_alloc = np.minimum(
                        resStorage_maxFracForIrrigation_CA * resStorageTotal_alloc,
                        demand_Segment / self.var.Water_conveyance_efficiency)  # [M3]

                    # fraction of water abstracted versus water available for total segment reservoir volumes
                    ResAbstractFactor = np.where(resStorageTotal_alloc > 0,
                                                 divideValues(act_bigLakeResAbst_alloc, resStorageTotal_alloc),
                                                 0)
                    # Compressed version needs to be corrected as above
                    ResAbstractFactorC = np.compress(self.var.compress_LR, ResAbstractFactor)
                    ResAbstractFactorC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C,
                                                     ResAbstractFactorC)

                    self.var.lakeStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
                    self.var.lakeVolumeM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC
                    self.var.lakeResStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
                    self.var.reservoirStorageM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC

                    self.var.abstractedLakeReservoirM3C += self.var.reservoirStorageM3C * ResAbstractFactorC
                    np.put(self.var.abstractedLakeReservoirM3, self.var.decompress_LR,
                           self.var.abstractedLakeReservoirM3C)

                    self.var.lakeResStorage = globals.inZero.copy()
                    np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)

                    metRemainSegment = np.where(demand_Segment > 0,
                                                divideValues(act_bigLakeResAbst_alloc * self.var.Water_conveyance_efficiency,
                                                             demand_Segment), 0)  # by definition <= 1

                    self.var.act_bigLakeResAbst += remainNeedPre * metRemainSegment
                    self.var.act_SurfaceWaterAbstract += remainNeedPre * metRemainSegment

                    self.var.act_ResAbst = remainNeedPre * metRemainSegment

                    self.var.Res_Domestic = np.minimum(self.var.act_ResAbst,
                                                       pot_Res_Domestic)
                    self.var.Res_Livestock = np.minimum(self.var.act_ResAbst - self.var.Res_Domestic,
                                                        pot_Res_Livestock)
                    self.var.Res_Industry = np.minimum(
                        self.var.act_ResAbst - self.var.Res_Domestic - self.var.Res_Livestock,
                        pot_Res_Industry)

                # If sector- and source-specific abstractions are activated, then domestic, industrial, and
                #  livestock demands were attempted to be satisfied in the previous step. Otherwise, total demands
                #  not satisfied by previous sources is attempted.
                #
                # The remaining demand within each command area [M3] is put into a map where each cell in the
                # command area holds this total demand
                demand_Segment = np.where(self.var.reservoir_command_areas > 0,
                                          npareatotal(remainNeed2 * self.var.cellArea,
                                                      self.var.reservoir_command_areas),
                                          0)  # [M3]

                ## Reservoir associated with the Command Area
                #
                # If there is more than one reservoir in a command area,
                #   the storage of the reservoir with maximum storage in this time-step is chosen.
                # The map resStorageTotal_alloc holds this maximum reservoir storage
                #   within a command area in all cells within that command area

                ReservoirsThatAreCurrentlyReservoirs = np.where(self.var.waterBodyTypCTemp == 2,
                                                                self.var.reservoirStorageM3C, 0)
                reservoirStorageM3 = globals.inZero.copy()
                # np.put(reservoirStorageM3, self.var.decompress_LR, self.var.reservoirStorageM3C)
                np.put(reservoirStorageM3, self.var.decompress_LR, ReservoirsThatAreCurrentlyReservoirs)

                resStorageTotal_alloc = np.where(self.var.reservoir_command_areas > 0,
                                                 npareamaximum(reservoirStorageM3,
                                                               self.var.reservoir_command_areas), 0)  # [M3]

                # In the map resStorageTotal_allocC, the maximum storage from each allocation segment
                #   is held in all reservoir cells within that allocation segment.
                # We now correct to remove the reservoirs
                #   that are not this maximum-storage-reservoir for the command area.
                resStorageTotal_allocC = np.compress(self.var.compress_LR, resStorageTotal_alloc)
                resStorageTotal_allocC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C,
                                                     resStorageTotal_allocC)

                # The rules for the maximum amount of water to be released for irrigation
                #   are found for the chosen maximum-storage reservoir in each command area
                day_of_year = globals.dateVar['currDate'].timetuple().tm_yday

                if 'Reservoir_releases' in binding:
                    # resStorage_maxFracForIrrigation = 0.5 + globals.inZero.copy()
                    resStorage_maxFracForIrrigation = readnetcdf2('Reservoir_releases', day_of_year,
                                                                  useDaily='DOY', value='Fraction of Volume')
                elif self.var.reservoir_releases_excel_option:
                    resStorage_maxFracForIrrigation = globals.inZero.copy()
                    resStorage_maxFracForIrrigationC = np.where(self.var.lakeResStorage_release_ratioC > -1,
                                                                self.var.reservoir_supply[dateVar['doy']-1],
                                                                0.03)
                    np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR, resStorage_maxFracForIrrigationC)
                else:
                    resStorage_maxFracForIrrigation = 0.03 + globals.inZero.copy()

                # resStorage_maxFracForIrrigationC holds the fractional rules found for each reservoir,
                #   so we must null those that are not the maximum-storage reservoirs
                resStorage_maxFracForIrrigationC = np.compress(self.var.compress_LR, resStorage_maxFracForIrrigation)
                resStorage_maxFracForIrrigationC = np.multiply(
                    resStorageTotal_allocC == self.var.reservoirStorageM3C, resStorage_maxFracForIrrigationC)
                np.put(resStorage_maxFracForIrrigation, self.var.decompress_LR, resStorage_maxFracForIrrigationC)

                resStorage_maxFracForIrrigation_CA = np.where(self.var.reservoir_command_areas > 0,
                                                              npareamaximum(resStorage_maxFracForIrrigation,
                                                                            self.var.reservoir_command_areas), 0)

                act_bigLakeResAbst_alloc = np.minimum(resStorage_maxFracForIrrigation_CA * resStorageTotal_alloc,
                                                      demand_Segment / self.var.Water_conveyance_efficiency)  # [M3]

                ResAbstractFactor = np.where(resStorageTotal_alloc > 0,
                                             divideValues(act_bigLakeResAbst_alloc, resStorageTotal_alloc),
                                             0)
                # fraction of water abstracted versus water available for total segment reservoir volumes
                # Compressed version needs to be corrected as above
                ResAbstractFactorC = np.compress(self.var.compress_LR, ResAbstractFactor)
                ResAbstractFactorC = np.multiply(resStorageTotal_allocC == self.var.reservoirStorageM3C,
                                                 ResAbstractFactorC)

                self.var.lakeStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
                self.var.lakeVolumeM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC
                self.var.lakeResStorageC -= self.var.reservoirStorageM3C * ResAbstractFactorC
                self.var.reservoirStorageM3C -= self.var.reservoirStorageM3C * ResAbstractFactorC

                self.var.abstractedLakeReservoirM3C += self.var.reservoirStorageM3C * ResAbstractFactorC
                np.put(self.var.abstractedLakeReservoirM3, self.var.decompress_LR,
                       self.var.abstractedLakeReservoirM3C)

                self.var.lakeResStorage = globals.inZero.copy()
                np.put(self.var.lakeResStorage, self.var.decompress_LR, self.var.lakeResStorageC)

                metRemainSegment = np.where(demand_Segment > 0,
                                            divideValues(act_bigLakeResAbst_alloc * self.var.Water_conveyance_efficiency,
                                                         demand_Segment), 0)  # by definition <= 1

                self.var.leakageC_daily = resStorageTotal_allocC * ResAbstractFactorC * (
                        1 - np.compress(self.var.compress_LR, self.var.Water_conveyance_efficiency))

                self.var.leakage = globals.inZero.copy()
                np.put(self.var.leakage, self.var.decompress_LR, self.var.leakageC_daily)

                # self.var.leakageC += self.var.leakageC_daily
                divleak_canal = divideValues((self.var.leakageC_daily) ,self.var.canalsAreaC)
                self.var.leakageCanalsC_M = np.where(self.var.canalsAreaC > 0,divleak_canal, 0)

                # Without this, npareamaximum uses the historical maximum
                self.var.leakageCanals_M = globals.inZero.copy()
                np.put(self.var.leakageCanals_M, self.var.decompress_LR, self.var.leakageCanalsC_M)  # good
                self.var.leakageCanals_M = npareamaximum(self.var.leakageCanals_M,
                                                         self.var.canals)

                self.var.act_bigLakeResAbst += remainNeed2 * metRemainSegment
                self.var.act_SurfaceWaterAbstract += remainNeed2 * metRemainSegment

                self.var.act_ResAbst += remainNeed2 * metRemainSegment

                ## End of using_reservoir_command_areas

                if self.var.sectorSourceAbstractionFractions:
                    self.var.Res_Irrigation = np.minimum(
                        remainNeed2 * metRemainSegment,
                        pot_Res_Irrigation)

                # B

            # remaining is taken from groundwater if possible
            if self.var.sectorSourceAbstractionFractions:
                pot_GW_Domestic = np.minimum(
                    self.var.gwAbstractionFraction_Domestic * self.var.domesticDemand,
                    self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic - \
                    self.var.Lift_Domestic - self.var.Lake_Domestic - self.var.Res_Domestic)
                    
                pot_GW_Livestock = np.minimum(
                    self.var.gwAbstractionFraction_Livestock * self.var.livestockDemand,
                    self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock - \
                    self.var.Lift_Livestock - self.var.Lake_Livestock - self.var.Res_Livestock)

                pot_GW_Industry = np.minimum(
                    self.var.gwAbstractionFraction_Industry * self.var.industryDemand,
                    self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry - \
                    self.var.Lift_Industry - self.var.Lake_Industry - self.var.Res_Industry)

                pot_GW_Irrigation = np.minimum(
                    self.var.gwAbstractionFraction_Irrigation * self.var.totalIrrDemand,                    
                    self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation - \
                    self.var.Lift_Irrigation - self.var.Lake_Irrigation - self.var.Res_Irrigation)


                self.var.pot_GroundwaterAbstract = pot_GW_Domestic + pot_GW_Livestock + pot_GW_Industry + pot_GW_Irrigation
            else:
                self.var.pot_GroundwaterAbstract = totalDemand - self.var.act_SurfaceWaterAbstract

            self.var.nonFossilGroundwaterAbs = np.maximum(0., np.minimum(self.var.readAvlStorGroundwater,
                                                                             self.var.pot_GroundwaterAbstract))

            # if limitAbstraction from groundwater is True
            # fossil gwAbstraction and water demand may be reduced
            # variable to reduce/limit groundwater abstraction (> 0 if limitAbstraction = True)

            if self.var.sectorSourceAbstractionFractions:
                # A
                self.var.GW_Domestic = np.minimum(self.var.nonFossilGroundwaterAbs, pot_GW_Domestic)
                self.var.GW_Livestock = np.minimum(self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic,
                                                   pot_GW_Livestock)
                self.var.GW_Industry = np.minimum(
                    self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic - self.var.GW_Livestock,
                    pot_GW_Industry)
                self.var.GW_Irrigation = np.minimum(
                    self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic - self.var.GW_Livestock - self.var.GW_Industry,
                    pot_GW_Irrigation)

                unmet_Domestic = self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic - \
                    self.var.Lift_Domestic - self.var.Lake_Domestic - self.var.Res_Domestic - self.var.GW_Domestic
                unmet_Livestock = self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock - \
                    self.var.Lift_Livestock - self.var.Lake_Livestock - self.var.Res_Livestock - self.var.GW_Livestock
                unmet_Industry = self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry - \
                    self.var.Lift_Industry  - self.var.Lake_Industry - self.var.Res_Industry - self.var.GW_Industry
                unmet_Irrigation = self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation - \
                    self.var.Lift_Industry - self.var.Lake_Irrigation - self.var.Res_Irrigation - self.var.GW_Irrigation

            if checkOption('limitAbstraction'):
                # real surface water abstraction can be lower, because not all demand can be done from surface water
                act_swAbstractionFraction = divideValues(self.var.act_SurfaceWaterAbstract, totalDemand)
                # Fossil groundwater abstraction is not allowed
                # allocation rule here: domestic& industry > irrigation > paddy

                if self.var.sectorSourceAbstractionFractions:

                    self.var.act_nonIrrWithdrawal = self.var.Desal_Domestic + self.var.Desal_Livestock + self.var.Desal_Industry + \
                                                    self.var.Channel_Domestic + self.var.Channel_Livestock + self.var.Channel_Industry + \
                                                    self.var.Lift_Domestic + self.var.Lift_Livestock + self.var.Lift_Industry + \
                                                    self.var.Lake_Domestic + self.var.Lake_Livestock + self.var.Lake_Industry + \
                                                    self.var.Res_Domestic + self.var.Res_Livestock + self.var.Res_Industry + \
                                                    self.var.GW_Domestic + self.var.GW_Livestock + self.var.GW_Industry
                    self.var.act_irrWithdrawal = self.var.Desal_Irrigation + self.var.Channel_Irrigation + \
                        self.var.Lift_Irrigation + self.var.Lake_Irrigation + self.var.Res_Irrigation + self.var.GW_Irrigation
                    # Currently desalination is accounted as surface water
                    act_irrWithdrawalSW = self.var.Desal_Irrigation + self.var.Channel_Irrigation + self.var.Lift_Irrigation + \
                        self.var.Lake_Irrigation + self.var.Res_Irrigation
                    act_irrWithdrawalGW = self.var.GW_Irrigation
                    self.var.act_irrNonpaddyWithdrawal = np.minimum(self.var.act_irrWithdrawal,
                                                                    self.var.fracVegCover[3] * self.var.irrDemand[3])
                    self.var.act_irrPaddyWithdrawal = self.var.act_irrWithdrawal - self.var.act_irrNonpaddyWithdrawal

                    act_gw = np.copy(self.var.nonFossilGroundwaterAbs)

                elif self.var.includeIndusDomesDemand:  # all demands are taken into account
                    # non-irrgated water demand: adjusted (and maybe increased) by gwabstration factor if
                    # non-irrgated water demand is higher than actual growndwater abstraction (what is needed and
                    # what is stored in gw)
                    act_nonIrrWithdrawalGW = self.var.nonIrrDemand * (1 - act_swAbstractionFraction)
                    act_nonIrrWithdrawalGW = np.where(act_nonIrrWithdrawalGW > self.var.nonFossilGroundwaterAbs,
                                                      self.var.nonFossilGroundwaterAbs, act_nonIrrWithdrawalGW)
                    act_nonIrrWithdrawalSW = act_swAbstractionFraction * self.var.nonIrrDemand
                    self.var.act_nonIrrWithdrawal = act_nonIrrWithdrawalSW + act_nonIrrWithdrawalGW

                    # irrigated water demand:
                    act_irrWithdrawalGW = self.var.totalIrrDemand * (1 - act_swAbstractionFraction)
                    act_irrWithdrawalGW = np.minimum(self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW,
                                                     act_irrWithdrawalGW)
                    act_irrWithdrawalSW = act_swAbstractionFraction * self.var.totalIrrDemand
                    self.var.act_irrWithdrawal = act_irrWithdrawalSW + act_irrWithdrawalGW
                    # (nonpaddy)
                    act_irrnonpaddyGW = self.var.fracVegCover[3] * (1 - act_swAbstractionFraction) * \
                                        self.var.irrDemand[3]
                    act_irrnonpaddyGW = np.minimum(self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW,
                                                   act_irrnonpaddyGW)
                    act_irrnonpaddySW = self.var.fracVegCover[3] * act_swAbstractionFraction * self.var.irrDemand[3]
                    self.var.act_irrNonpaddyWithdrawal = act_irrnonpaddySW + act_irrnonpaddyGW
                    # (paddy)
                    act_irrpaddyGW = self.var.fracVegCover[2] * (1 - act_swAbstractionFraction) * self.var.irrDemand[2]
                    act_irrpaddyGW = np.minimum(
                        self.var.nonFossilGroundwaterAbs - act_nonIrrWithdrawalGW - act_irrnonpaddyGW, act_irrpaddyGW)
                    act_irrpaddySW = self.var.fracVegCover[2] * act_swAbstractionFraction * self.var.irrDemand[2]
                    self.var.act_irrPaddyWithdrawal = act_irrpaddySW + act_irrpaddyGW

                    act_gw = act_nonIrrWithdrawalGW + act_irrWithdrawalGW
                    # This should be equal to self.var.nonFossilGroundwaterAbs?


                else:  # only irrigation is considered

                    self.var.act_nonIrrWithdrawal = globals.inZero.copy()

                    # irrigated water demand:
                    act_irrWithdrawalGW = self.var.totalIrrDemand * (1 - act_swAbstractionFraction)
                    act_irrWithdrawalGW = np.minimum(self.var.nonFossilGroundwaterAbs, act_irrWithdrawalGW)
                    act_irrWithdrawalSW = act_swAbstractionFraction * self.var.totalIrrDemand
                    self.var.act_irrWithdrawal = act_irrWithdrawalSW + act_irrWithdrawalGW
                    # (nonpaddy)
                    act_irrnonpaddyGW = self.var.fracVegCover[3] * (1 - act_swAbstractionFraction) * \
                                        self.var.irrDemand[3]
                    act_irrnonpaddyGW = np.minimum(self.var.nonFossilGroundwaterAbs, act_irrnonpaddyGW)
                    act_irrnonpaddySW = self.var.fracVegCover[3] * act_swAbstractionFraction * self.var.irrDemand[3]
                    self.var.act_irrNonpaddyWithdrawal = act_irrnonpaddySW + act_irrnonpaddyGW
                    # (paddy)
                    act_irrpaddyGW = self.var.fracVegCover[2] * (1 - act_swAbstractionFraction) * self.var.irrDemand[2]
                    act_irrpaddyGW = np.minimum(
                        self.var.nonFossilGroundwaterAbs - act_irrnonpaddyGW, act_irrpaddyGW)
                    act_irrpaddySW = self.var.fracVegCover[2] * act_swAbstractionFraction * self.var.irrDemand[2]
                    self.var.act_irrPaddyWithdrawal = act_irrpaddySW + act_irrpaddyGW

                    act_gw = np.copy(act_irrWithdrawalGW)

                # calculate act_ water demand, because irr demand has still demand from previous day included
                # if the demand from previous day is not fulfilled it is taken to the next day and so on
                # if we do not correct we double account each day the demand from previous days
                self.var.act_irrPaddyDemand = np.maximum(0, self.var.irrPaddyDemand - self.var.unmetDemandPaddy)
                self.var.act_irrNonpaddyDemand = np.maximum(0,
                                                            self.var.irrNonpaddyDemand - self.var.unmetDemandNonpaddy)

                # unmet is either pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs or demand - withdrawal
                if self.var.includeIndusDomesDemand:  # all demands are taken into account
                    self.var.unmetDemand = (self.var.totalIrrDemand - self.var.act_irrWithdrawal) + \
                                           (self.var.nonIrrDemand - self.var.act_nonIrrWithdrawal)
                else:  # only irrigation is considered
                    self.var.unmetDemand = (self.var.totalIrrDemand - self.var.act_irrWithdrawal) - \
                                           self.var.act_nonIrrWithdrawal
                self.var.unmetDemandPaddy = self.var.irrPaddyDemand - self.var.act_irrPaddyDemand
                self.var.unmetDemandNonpaddy = self.var.irrNonpaddyDemand - self.var.act_irrNonpaddyDemand


            else:
                # Fossil groundwater abstractions are allowed (act = pot)
                if 'zonal_abstraction' in option:
                    if checkOption('zonal_abstraction'):

                        # using allocation from abstraction zone
                        # this might be a regular grid e.g. 2x2 for 0.5 deg
                        left_sf = self.var.readAvlChannelStorageM  # already removed - self.var.act_channelAbst
                        # sum demand, surface water - local used, groundwater - local use, not satisfied for allocation zone

                        if self.var.sectorSourceAbstractionFractions:
                            unmetChannel_Domestic = pot_Channel_Domestic - self.var.Channel_Domestic
                            unmetChannel_Livestock = pot_Channel_Livestock - self.var.Channel_Livestock
                            unmetChannel_Industry = pot_Channel_Industry - self.var.Channel_Industry
                            unmetChannel_Irrigation = pot_Channel_Irrigation - self.var.Channel_Irrigation

                            pot_Channel_Domestic = np.minimum(unmetChannel_Domestic, unmet_Domestic)
                            pot_Channel_Livestock = np.minimum(unmetChannel_Livestock, unmet_Livestock)
                            pot_Channel_Industry = np.minimum(unmetChannel_Industry, unmet_Industry)
                            pot_Channel_Irrigation = np.minimum(unmetChannel_Irrigation, unmet_Irrigation)

                            unmet_Channel = pot_Channel_Domestic + pot_Channel_Livestock \
                                            + pot_Channel_Industry + pot_Channel_Irrigation

                            zoneDemand = npareatotal(unmet_Channel * self.var.cellArea, self.var.allocation_zone)

                        else:
                            zoneDemand = npareatotal(self.var.unmetDemand * self.var.cellArea, self.var.allocation_zone)

                        zone_sf_avail = npareatotal(left_sf, self.var.allocation_zone)

                        # zone abstraction is minimum of availability and demand
                        zone_sf_abstraction = np.minimum(zoneDemand, zone_sf_avail)
                        # water taken from surface zone and allocated to cell demand
                        cell_sf_abstraction = np.maximum(0., divideValues(left_sf, zone_sf_avail) * zone_sf_abstraction)
                        cell_sf_allocation = np.maximum(0., divideValues(self.var.unmetDemand,
                                                                            zoneDemand) * zone_sf_abstraction)

                        # sum up with other abstraction
                        self.var.act_SurfaceWaterAbstract = self.var.act_SurfaceWaterAbstract + cell_sf_abstraction
                        self.var.act_channelAbst = self.var.act_channelAbst + cell_sf_abstraction

                        if self.var.sectorSourceAbstractionFractions:
                            self.var.Channel_Domestic_fromZone = np.minimum(cell_sf_abstraction, pot_Channel_Domestic)
                            self.var.Channel_Domestic += self.var.Channel_Domestic_fromZone

                            self.var.Channel_Livestock_fromZone = np.minimum(
                                cell_sf_abstraction - self.var.Channel_Domestic_fromZone,
                                pot_Channel_Livestock)
                            self.var.Channel_Livestock += self.var.Channel_Livestock_fromZone

                            self.var.Channel_Industry_fromZone = np.minimum(
                                cell_sf_abstraction - self.var.Channel_Domestic_fromZone -
                                self.var.Channel_Livestock_fromZone,
                                pot_Channel_Industry)
                            self.var.Channel_Industry += self.var.Channel_Industry_fromZone

                            self.var.Channel_Irrigation_fromZone = np.minimum(
                                cell_sf_abstraction - self.var.Channel_Domestic_fromZone -
                                self.var.Channel_Livestock_fromZone - self.var.Channel_Industry_fromZone,
                                pot_Channel_Irrigation)
                            self.var.Channel_Irrigation += self.var.Channel_Irrigation_fromZone

                        # new potential groundwater abstraction
                        self.var.pot_GroundwaterAbstract = \
                            np.maximum(0., self.var.pot_GroundwaterAbstract - cell_sf_allocation)

                        left_gw_demand = np.maximum(0., self.var.pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs)
                        left_gw_avail = self.var.readAvlStorGroundwater - self.var.nonFossilGroundwaterAbs
                        zone_gw_avail = npareatotal(left_gw_avail * self.var.cellArea, self.var.allocation_zone)

                        # for groundwater substract demand which is fulfilled by surface zone, calc abstraction and what
                        # is left. zone_gw_demand = npareatotal(left_gw_demand, self.var.allocation_zone)
                        zone_gw_demand = zoneDemand - zone_sf_abstraction
                        zone_gw_abstraction = np.minimum(zone_gw_demand, zone_gw_avail)
                        # zone_unmetdemand = np.maximum(0., zone_gw_demand - zone_gw_abstraction)

                        # water taken from groundwater zone and allocated to cell demand
                        cell_gw_abstraction = \
                            np.maximum(0., divideValues(left_gw_avail, zone_gw_avail) * zone_gw_abstraction)
                        cell_gw_allocation = \
                            np.maximum(0., divideValues(left_gw_demand, zone_gw_demand) * zone_gw_abstraction)

                        self.var.unmetDemand = np.maximum(0., left_gw_demand - cell_gw_allocation)
                        self.var.nonFossilGroundwaterAbs = self.var.nonFossilGroundwaterAbs + cell_gw_abstraction

                        # UNDER CONSTRUCTION
                        if self.var.sectorSourceAbstractionFractions:
                            self.var.GW_Domestic_fromZone = np.minimum(self.var.nonFossilGroundwaterAbs, pot_GW_Domestic)
                            self.var.GW_Domestic += self.var.GW_Domestic_fromZone.copy()

                            self.var.GW_Livestock_fromZone = np.minimum(
                                self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic_fromZone,
                                pot_GW_Livestock)
                            self.var.GW_Livestock += self.var.GW_Livestock_fromZone.copy()

                            self.var.GW_Industry_fromZone = np.minimum(
                                self.var.nonFossilGroundwaterAbs -
                                self.var.GW_Domestic_fromZone - self.var.GW_Livestock_fromZone,
                                pot_GW_Industry)
                            self.var.GW_Industry += self.var.GW_Industry_fromZone.copy()

                            self.var.GW_Irrigation_fromZone = np.minimum(
                                self.var.nonFossilGroundwaterAbs - self.var.GW_Domestic_fromZone -
                                self.var.GW_Livestock_fromZone - self.var.GW_Industry_fromZone,
                                pot_GW_Irrigation)
                            self.var.GW_Irrigation += self.var.GW_Irrigation_fromZone.copy()

                        # end of zonal abstraction

                    else:
                        self.var.unmetDemand = self.var.pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs
                        if self.var.sectorSourceAbstractionFractions:
                            self.var.GW_Domestic = pot_GW_Domestic
                            self.var.GW_Industry = pot_GW_Industry
                            self.var.GW_Livestock = pot_GW_Livestock
                            self.var.GW_Irrigation = pot_GW_Irrigation

                else:
                    self.var.unmetDemand = self.var.pot_GroundwaterAbstract - self.var.nonFossilGroundwaterAbs
                    if self.var.sectorSourceAbstractionFractions:
                        self.var.GW_Domestic = pot_GW_Domestic
                        self.var.GW_Industry = pot_GW_Industry
                        self.var.GW_Livestock = pot_GW_Livestock
                        self.var.GW_Irrigation = pot_GW_Irrigation


                if self.var.includeIndusDomesDemand:  # all demands are taken into account
                    self.var.act_nonIrrWithdrawal = np.copy(self.var.nonIrrDemand)
                self.var.act_irrWithdrawal = np.copy(self.var.totalIrrDemand)

                act_gw = np.copy(self.var.pot_GroundwaterAbstract)

                self.var.act_irrNonpaddyWithdrawal = self.var.fracVegCover[3] * self.var.irrDemand[3]
                self.var.act_irrPaddyWithdrawal = self.var.fracVegCover[2] * self.var.irrDemand[2]

                self.var.act_irrNonpaddyDemand = self.var.act_irrNonpaddyWithdrawal.copy()
                self.var.act_irrPaddyDemand = self.var.act_irrPaddyWithdrawal.copy()

            ## End of limit extraction if, then

            self.var.act_irrConsumption[2] = divideValues(self.var.act_irrPaddyWithdrawal,
                                                          self.var.fracVegCover[2]) * self.var.efficiencyPaddy
            self.var.act_irrConsumption[3] = divideValues(self.var.act_irrNonpaddyWithdrawal,
                                                          self.var.fracVegCover[3]) * self.var.efficiencyNonpaddy

            if self.var.sectorSourceAbstractionFractions:

                self.var.act_domWithdrawal = self.var.Channel_Domestic + self.var.Lift_Domestic + \
                                             self.var.Desal_Domestic + self.var.Lake_Domestic + \
                                             self.var.Res_Domestic + self.var.GW_Domestic
                self.var.act_livWithdrawal = self.var.Channel_Livestock + self.var.Lift_Livestock + \
                                             self.var.Desal_Livestock + self.var.Lake_Livestock + \
                                             self.var.Res_Livestock + self.var.GW_Livestock
                self.var.act_indWithdrawal = self.var.Channel_Industry + self.var.Lift_Industry + \
                                             self.var.Desal_Industry + self.var.Lake_Industry + \
                                             self.var.Res_Industry + self.var.GW_Industry

                self.var.act_indConsumption = self.var.ind_efficiency * self.var.act_indWithdrawal
                self.var.act_domConsumption = self.var.dom_efficiency * self.var.act_domWithdrawal
                self.var.act_livConsumption = self.var.liv_efficiency * self.var.act_livWithdrawal
                self.var.act_nonIrrConsumption = self.var.act_domConsumption + self.var.act_indConsumption + \
                                                 self.var.act_livConsumption

            elif self.var.includeIndusDomesDemand:  # all demands are taken into account

                self.var.act_indWithdrawal = frac_industry * self.var.act_nonIrrWithdrawal
                self.var.act_domWithdrawal = frac_domestic * self.var.act_nonIrrWithdrawal
                self.var.act_livWithdrawal = frac_livestock * self.var.act_nonIrrWithdrawal
                self.var.act_indConsumption = self.var.ind_efficiency * self.var.act_indWithdrawal
                self.var.act_domConsumption = self.var.dom_efficiency * self.var.act_domWithdrawal
                self.var.act_livConsumption = self.var.liv_efficiency * self.var.act_livWithdrawal
                self.var.act_nonIrrConsumption = self.var.act_domConsumption + self.var.act_indConsumption + \
                                                 self.var.act_livConsumption

            else:  # only irrigation is considered
                self.var.act_nonIrrConsumption = globals.inZero.copy()

            self.var.act_totalIrrConsumption = self.var.fracVegCover[2] * self.var.act_irrConsumption[2] + \
                                               self.var.fracVegCover[3] * self.var.act_irrConsumption[3]
            self.var.act_paddyConsumption = self.var.fracVegCover[2] * self.var.act_irrConsumption[2]
            self.var.act_nonpaddyConsumption = self.var.fracVegCover[3] * self.var.act_irrConsumption[3]

            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.totalWaterDemand = self.var.fracVegCover[2] * self.var.irrDemand[2] + self.var.fracVegCover[
                    3] * self.var.irrDemand[3] + self.var.nonIrrDemand
                self.var.act_totalWaterWithdrawal = self.var.act_nonIrrWithdrawal + self.var.act_irrWithdrawal
                self.var.act_totalWaterConsumption = self.var.act_nonIrrConsumption + self.var.act_totalIrrConsumption
            else:  # only irrigation is considered
                self.var.totalWaterDemand = self.var.fracVegCover[2] * self.var.irrDemand[2] + self.var.fracVegCover[
                    3] * self.var.irrDemand[3]
                self.var.act_totalWaterWithdrawal = np.copy(self.var.act_irrWithdrawal)
                self.var.act_totalWaterConsumption = np.copy(self.var.act_totalIrrConsumption)

            # --- calculate return flow
            # Sum up loss - difference between withdrawn and consumed - split into return flow and evaporation
            sumIrrLoss = self.var.act_irrWithdrawal - self.var.act_totalIrrConsumption

            self.var.returnflowIrr = self.var.returnfractionIrr * sumIrrLoss
            self.var.addtoevapotrans = (1 - self.var.returnfractionIrr) * sumIrrLoss

            if self.var.sectorSourceAbstractionFractions:
                self.var.returnflowNonIrr = self.var.act_nonIrrWithdrawal - self.var.act_nonIrrConsumption
            elif self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.returnflowNonIrr = self.var.nonIrrReturnFlowFraction * self.var.act_nonIrrWithdrawal

            # limit return flow to not put all fossil groundwater back into the system, because it can lead to higher
            # river discharge than without water demand, as water is taken from fossil groundwater (out of system)
            unmet_div_ww = 1. - np.minimum(1, divideValues(self.var.unmetDemand,
                                                           self.var.act_totalWaterWithdrawal + self.var.unmetDemand))

            # 'fossil_water_treated_normally' means that there is no lost fossil water
            if 'fossil_water_treated_normally' in option:
                if checkOption('fossil_water_treated_normally'):
                    unmet_div_ww = 1

            if checkOption('limitAbstraction'):
                unmet_div_ww = 1

            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.unmet_lost = (self.var.returnflowIrr + self.var.returnflowNonIrr + self.var.addtoevapotrans) \
                                      * (1 - unmet_div_ww)
            else:  # only irrigation is considered
                self.var.unmet_lost = (self.var.returnflowIrr + self.var.addtoevapotrans) * (1 - unmet_div_ww)

            # self.var.waterDemandLost = self.var.act_totalWaterConsumption + self.var.addtoevapotrans
            self.var.unmet_lostirr = (self.var.returnflowIrr + self.var.addtoevapotrans) * (1 - unmet_div_ww)
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.unmet_lostNonirr = self.var.returnflowNonIrr * (1 - unmet_div_ww)

            self.var.returnflowIrr = self.var.returnflowIrr * unmet_div_ww
            self.var.addtoevapotrans = self.var.addtoevapotrans * unmet_div_ww
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.returnflowNonIrr = self.var.returnflowNonIrr * unmet_div_ww            

            # returnflow to river and to evapotranspiration
            if self.var.includeIndusDomesDemand:  # all demands are taken into account
                self.var.returnFlow = self.var.returnflowIrr + self.var.returnflowNonIrr
            else:  # only irrigation is considered
                self.var.returnFlow = self.var.returnflowIrr

            self.var.waterabstraction = self.var.nonFossilGroundwaterAbs + self.var.unmetDemand + \
                                        self.var.act_SurfaceWaterAbstract


    # ---------------------------------------------------------------------------
    # ===========================================================================
    # ---------------------------------------------------------------------------

    def calc_totalDemand(self):    
        """
        Calculate total water demand (irrigation + non-irrigation) and fractions of non-irrigation demand 
        """
        if self.var.includeIndusDomesDemand:  # all demands are taken into account
            if globals.dateVar['newStart'] or globals.dateVar['newMonth'] \
                    or 'reservoir_transfers' in option:
                # total (potential) non irrigation water demand
                self.var.nonIrrDemand = self.var.domesticDemand + self.var.industryDemand + self.var.livestockDemand
                self.var.pot_nonIrrConsumption = np.minimum(self.var.nonIrrDemand,
                                                            self.var.pot_domesticConsumption +
                                                            self.var.pot_industryConsumption + self.var.pot_livestockConsumption)
                # fraction of return flow from domestic and industrial water demand
                self.var.nonIrrReturnFlowFraction = divideValues(
                    (self.var.nonIrrDemand - self.var.pot_nonIrrConsumption), self.var.nonIrrDemand)

            # non-irrg fracs in nonIrrDemand
            frac_industry = divideValues(self.var.industryDemand, self.var.nonIrrDemand)
            frac_domestic = divideValues(self.var.domesticDemand, self.var.nonIrrDemand)
            frac_livestock = divideValues(self.var.livestockDemand, self.var.nonIrrDemand)

            # Sum up water demand
            # totalDemand [m]: total maximum (potential) water demand: irrigation and non irrigation
            totalDemand = self.var.nonIrrDemand + self.var.totalIrrDemand  # in [m]
        else:  # only irrigation is considered
            totalDemand = np.copy(self.var.totalIrrDemand)  # in [m]
            self.var.nonIrrDemand = globals.inZero.copy()
            self.var.pot_nonIrrConsumption = globals.inZero.copy()
            self.var.nonIrrReturnFlowFraction = globals.inZero.copy()
            frac_industry = globals.inZero.copy()
            frac_domestic = globals.inZero.copy()
            frac_livestock = globals.inZero.copy()

            return totalDemand, frac_industry, frac_domestic, frac_livestock   


    def calc_desalination(self):
        """
        Calculate desalination abstraction and allocation to sectors
        """  
        self.var.act_DesalWaterAbstractM3 = globals.inZero.copy()
        # Desalination is not allowed without sectorSourceAbstractionFractions
        if self.var.sectorSourceAbstractionFractions:
            if self.var.includeDesal:
                pot_Desal_Domestic = self.var.othAbstractionFraction_Desal_Domestic * self.var.domesticDemand
                pot_Desal_Livestock = self.var.othAbstractionFraction_Desal_Livestock * self.var.livestockDemand
                pot_Desal_Industry = self.var.othAbstractionFraction_Desal_Industry * self.var.industryDemand
                pot_Desal_Irrigation = self.var.othAbstractionFraction_Desal_Irrigation * self.var.totalIrrDemand

                pot_DesalAbst = pot_Desal_Domestic + pot_Desal_Livestock + pot_Desal_Industry + pot_Desal_Irrigation
                if not self.var.unlimitedDesal:
                    self.var.AvlDesalM3 = self.var.desalAnnualCap[dateVar['currDate'].year] / 365
                    abstractLimitCoeff = np.minimum(np.nansum(pot_DesalAbst * self.var.cellArea), self.var.AvlDesalM3) / np.nansum(pot_DesalAbst * self.var.cellArea)
                    self.var.act_DesalWaterAbstractM = pot_DesalAbst * abstractLimitCoeff
                else:
                    self.var.act_DesalWaterAbstractM = pot_DesalAbst
    
                self.var.Desal_Domestic = np.minimum(self.var.act_DesalWaterAbstractM,
                                                    self.var.othAbstractionFraction_Desal_Domestic * self.var.domesticDemand)
                self.var.Desal_Livestock = np.minimum(self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic,
                                                        self.var.othAbstractionFraction_Desal_Livestock * self.var.livestockDemand)
                self.var.Desal_Industry = np.minimum(
                    self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic - self.var.Desal_Livestock,
                    self.var.othAbstractionFraction_Desal_Industry * self.var.industryDemand)
                self.var.Desal_Irrigation = np.minimum(
                    self.var.act_DesalWaterAbstractM - self.var.Desal_Domestic - self.var.Desal_Livestock - self.var.Desal_Industry,
                    self.var.othAbstractionFraction_Desal_Irrigation * self.var.totalIrrDemand)


    def calc_channelAbstract_Lift(self):
        """
        Calculate actual channel abstraction from lift areas
        """
        # UNDER CONSTRUCTION
        if self.var.using_lift_areas:
            # Lift development When there is sufficient water in the Segment to fulfill demand, the water is
            # taken away proportionally from each cell's readAvlChannelStorageM in the Segment. For example,
            # if total demand can be filled with 50% of total availability, then 50% of the
            # readAvlChannelStorageM from each cell is used. Note that if a cell has too little Channel Storage,
            # then no water will be taken from the cell as this was dealt with earlier:  readAvlChannelStorage =
            # 0 if < (0.0005 * self.var.cellArea) Note: Due to the shared use of abstracted channel storage,
            # a cell may abstract more than its pot_SurfaceAbstract, as well as not necessarily satisfy its
            # pot_SurfaceAbstract
            
            pot_Lift_Domestic = np.minimum(self.var.swAbstractionFraction_Lift_Domestic * self.var.domesticDemand, \
                self.var.domesticDemand - self.var.Desal_Domestic - self.var.Channel_Domestic )                
            pot_Lift_Livestock = np.minimum(self.var.swAbstractionFraction_Lift_Livestock * self.var.livestockDemand, \
                self.var.livestockDemand - self.var.Desal_Livestock - self.var.Channel_Livestock)
            pot_Lift_Industry = np.minimum(self.var.swAbstractionFraction_Lift_Industry * self.var.industryDemand, \
                self.var.industryDemand - self.var.Desal_Industry - self.var.Channel_Industry )
            pot_Lift_Irrigation = np.minimum(self.var.swAbstractionFraction_Lift_Irrigation * self.var.totalIrrDemand, \
                self.var.totalIrrDemand - self.var.Desal_Irrigation - self.var.Channel_Irrigation)

            pot_liftAbst = pot_Lift_Domestic + pot_Lift_Livestock + pot_Lift_Industry + pot_Lift_Irrigation


            remainNeed_afterLocal = pot_liftAbst.copy()

            # The remaining demand within each command area [M3] is put into a map where each cell in the command
            # area holds this total demand
            demand_Segment_lift = np.where(self.var.lift_command_areas > 0,
                                            npareatotal(remainNeed_afterLocal * self.var.cellArea,
                                                        self.var.lift_command_areas),
                                            0)  # [M3]

            available_Segment_lift = np.where(self.var.lift_command_areas > 0,
                                                npareatotal(self.var.readAvlChannelStorageM * self.var.cellArea,
                                                            self.var.lift_command_areas),
                                                0)  # [M3]

            frac_used_Segment_lift = np.where(available_Segment_lift > 0,
                                                np.minimum(demand_Segment_lift / available_Segment_lift, 1.), 0.)

            self.var.act_channelAbst += (frac_used_Segment_lift * self.var.readAvlChannelStorageM)

            metRemainSegment_lift = np.where(demand_Segment_lift > 0,
                                                divideValues(frac_used_Segment_lift * available_Segment_lift,
                                                            demand_Segment_lift), 0)
            self.var.metRemainSegment_lift = metRemainSegment_lift.copy()
            lift_abstractions = metRemainSegment_lift * remainNeed_afterLocal
            self.var.act_SurfaceWaterAbstract += lift_abstractions
            self.var.readAvlChannelStorageM -= (frac_used_Segment_lift * self.var.readAvlChannelStorageM)
            self.var.readAvlChannelStorageM = np.where(self.var.readAvlChannelStorageM < 0.02, 0,
                                                        self.var.readAvlChannelStorageM)
            # Used in landCover for riverbedExchange
            self.var.act_channelAbstract_Lift = frac_used_Segment_lift * self.var.readAvlChannelStorageM

            if self.var.sectorSourceAbstractionFractions:

                # A
                self.var.Lift_Domestic = np.minimum(lift_abstractions, pot_Lift_Domestic)
                self.var.Lift_Livestock = np.minimum(lift_abstractions - self.var.Lift_Domestic,
                                                        pot_Lift_Livestock)
                self.var.Lift_Industry = np.minimum(
                    lift_abstractions - self.var.Lift_Domestic - self.var.Lift_Livestock,
                    pot_Lift_Industry)
                self.var.Lift_Irrigation = np.minimum(
                    lift_abstractions - self.var.Lift_Domestic - self.var.Lift_Livestock - self.var.Lift_Industry,
                    pot_Lift_Irrigation)

    def calc_channelAbstract_SourceAbstraction(self):
        """
        Calculate channel abstraction for given sector-specific source abstraction fractions 
        """
        if self.var.sectorSourceAbstractionFractions:                            
            pot_Channel_Domestic = np.minimum(self.var.swAbstractionFraction_Channel_Domestic * self.var.domesticDemand, \
                self.var.domesticDemand - self.var.Desal_Domestic)                
            pot_Channel_Livestock = np.minimum(self.var.swAbstractionFraction_Channel_Livestock * self.var.livestockDemand, \
                self.var.livestockDemand - self.var.Desal_Livestock)
            pot_Channel_Industry = np.minimum(self.var.swAbstractionFraction_Channel_Industry * self.var.industryDemand, \
                self.var.industryDemand - self.var.Desal_Industry)
            pot_Channel_Irrigation = np.minimum(self.var.swAbstractionFraction_Channel_Irrigation * self.var.totalIrrDemand, \
                self.var.totalIrrDemand - self.var.Desal_Irrigation)

            pot_channelAbst = pot_Channel_Domestic + pot_Channel_Livestock + pot_Channel_Industry + pot_Channel_Irrigation

            self.var.act_SurfaceWaterAbstract = np.minimum(self.var.readAvlChannelStorageM, pot_channelAbst)
        else:
            pot_SurfaceAbstract = totalDemand * self.var.swAbstractionFraction
            # only local surface water abstraction is allowed (network is only within a cell)
            self.var.act_SurfaceWaterAbstract = np.minimum(self.var.readAvlChannelStorageM, pot_SurfaceAbstract)

        self.var.readAvlChannelStorageM -= self.var.act_SurfaceWaterAbstract
        self.var.act_channelAbst = self.var.act_SurfaceWaterAbstract.copy()
        # if surface water is not sufficient it is taken from groundwater

        if self.var.sectorSourceAbstractionFractions:
            self.var.Channel_Domestic = np.minimum(self.var.act_channelAbst,
                                                    self.var.swAbstractionFraction_Channel_Domestic * self.var.domesticDemand)
            self.var.Channel_Livestock = np.minimum(self.var.act_channelAbst - self.var.Channel_Domestic,
                                                    self.var.swAbstractionFraction_Channel_Livestock * self.var.livestockDemand)
            self.var.Channel_Industry = np.minimum(
                self.var.act_channelAbst - self.var.Channel_Domestic - self.var.Channel_Livestock,
                self.var.swAbstractionFraction_Channel_Industry * self.var.industryDemand)
            self.var.Channel_Irrigation = np.minimum(
                self.var.act_channelAbst - self.var.Channel_Domestic - self.var.Channel_Livestock - self.var.Channel_Industry,
                self.var.swAbstractionFraction_Channel_Irrigation * self.var.totalIrrDemand)
