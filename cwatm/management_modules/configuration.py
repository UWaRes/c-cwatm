# -------------------------------------------------------------------------
# Name:        Configuration
# Purpose:
#
# Author:      burekpe
#
# Created:     16/05/2016
# Copyright:   (c) burekpe 2016
# -------------------------------------------------------------------------

from cwatm.management_modules.globals import *

import configparser
import re
import xml.dom.minidom
from cwatm.management_modules.messages import *

import os


import difflib  # to check the closest word in settingsfile, if an error occurs


class ExtParser(configparser.ConfigParser):
    """
    Extended configuration parser that supports placeholder substitution across sections.

    This class extends Python's built-in `ConfigParser` to allow dynamic referencing of values
    from other sections or within the same section using placeholder syntax.

    Supported placeholder formats:
    - `$(Section:Option)` — references an option from another section.
    - `$(Option)` — references an option from the same section.

    Example
    --------
        [FILE_PATHS]
        PathRoot = C:/work

        [INPUT]
        MaskMap = $(FILE_PATHS:PathRoot)/data/areamaps/area.tif

    Attributes
    ----------
    cur_depth : int
        Tracks the current depth of recursive interpolation to prevent infinite loops.

    as/copilot    
    """

    #implementing extended interpolation
    def __init__(self, *args, **kwargs):
        """
        Initializes the extended parser and sets the interpolation depth counter.
        """
        self.cur_depth = 0
        configparser.ConfigParser.__init__(self, *args, **kwargs)

    def get(self, section, option, raw=False, vars=None, **kwargs):
        """
        Retrieves the value of an option from a given section of the settingsfile, with support for 
        placeholder substitution.

        Parameters
        ----------
        section : str
            The section in the settings file.
        option : str
            The option (key) within the section.
        raw : bool, optional
            If True, disables interpolation and returns the raw value. Default is False.
        vars : dict, optional
            Additional variables for interpolation.

        Returns
        -------
        str
            The fully resolved value of the configuration option.

        Raises
        ------
        CWATMError
            If the requested option is not found and no close match is available.
        configparser.InterpolationDepthError
            If the maximum interpolation depth is exceeded due to recursive references.

        as/copilot
        """

        try:
           r_opt = configparser.ConfigParser.get(self, section, option, raw=True, vars=vars)
        except:
             print(section, option)
             closest = difflib.get_close_matches(option, list(binding.keys()))
             if not closest: closest = ["- no match -"]
             msg = "Error 116: Closest key to the required one is: \"" + closest[0] + "\""
             raise CWATMError(msg)

        if raw:
            return r_opt

        ret = r_opt
        re_newintp1 = r'\$\((\w*):(\w*)\)'  # other section
        re_newintp2 = r'\$\((\w*)\)'  # same section

        re_old1 = re.findall('\$\(\w*:\w*\)', r_opt)
        re_old2 = re.findall('\$\(\w*\)', r_opt)

        m_new1 = re.findall(re_newintp1, r_opt)
        m_new2 = re.findall(re_newintp2, r_opt)

        if m_new1:
             i = 0
             for f_section, f_option in m_new1:
                 self.cur_depth = self.cur_depth + 1
                 if self.cur_depth < configparser.MAX_INTERPOLATION_DEPTH:
                     sub = self.get(f_section,f_option, vars=vars)
                     ret = ret.replace(re_old1[i], sub)
                     i += 1
                 else:
                     raise configparser.InterpolationDepthError(option, section, r_opt)

        if m_new2:
             i = 0
             for l_option in m_new2:
                 self.cur_depth = self.cur_depth + 1
                 if self.cur_depth < configparser.MAX_INTERPOLATION_DEPTH:
                     sub = self.get(section, l_option, vars=vars)
                     ret = ret.replace(re_old2[i], sub)
                     i =+ 1
                 else:
                     raise configparser.InterpolationDepthError(option, section, r_opt)

        self.cur_depth = self.cur_depth - 1
        return ret



def parse_configuration(settingsFileName):
    """
    Parses the settings file and extracts simulation parameters, options, and output settings.

    This function reads a structured settings file using an extended parser (`ExtParser`) and organizes its contents
    into several global dictionaries and lists.

    Parameters
    ----------
    settingsFileName : str
        Path to the settings file to be parsed.

    Returns
    -------
    None
        This function modifies global variables in place and does not return a value.
        - binding: maps configuration keys to their values for simulation parameters.
        - option: stores boolean or integer flags from the OPTIONS section.
        - outMap and outTss: store output map and time series settings.
        - outDir: stores output directory paths.
        - outsection: tracks sections that contain output definitions.
        - outputDir: stores the main output directory path.

    Raises
    ------
    CWATMFileError
        If the specified settings file does not exist or cannot be found.

    Used in
    ------
        - run_cwatm.py (CWATMexe, CWATMexe2)
        
    as/copilot
    """

    def splitout(varin, check):
        """
        split variable in several one, seperator = ,

        :param varin:
        :param check:
        :return: list with several variables
        """

        out = list(map(str.strip, varin.split(',')))
        if out[0] == "": out[0]="None"
        if out[0] != "None": check = True
        return out, check

    if not(os.path.isfile(settingsFileName)):
        msg = "Error 302: Settingsfile not found!\n"
        raise CWATMFileError(settingsFileName,msg)
    config = ExtParser()
    config.optionxform = str
    config.sections()
    config.read(settingsFileName)
    for sec in config.sections():
        #print sec
        options = config.options(sec)
        check_section = False
        for opt in options:
            if sec=="OPTIONS":
                try:
                    option[opt] = config.getboolean(sec, opt)
                except:
                    option[opt] = config.getint(sec, opt)
            else:
                # Check if config line = output line
                if opt.lower()[0:4] == "out_":
                    index = sec.lower()+"_"+opt.lower()

                    if opt.lower()[-4:] =="_dir":
                        outDir[sec] = config.get(sec, opt)
                    else:
                        # split into timeseries and maps
                        if opt.lower()[4:8] == "tss_":
                            outTss[index],check_section = splitout(config.get(sec, opt),check_section)
                        else:
                            outMap[index],check_section = splitout(config.get(sec, opt),check_section)

                else:
                    # binding: all the parameters which are not output or option are collected
                    binding[opt] = config.get(sec, opt)

        if check_section:
            outsection.append(sec)

    outputDir.append(binding["PathOut"])
     # Output directory is stored in a separat global array


def read_metanetcdf(metaxml, name):
    """
    Read the metadata for netcdf output files
    unit, long name, standard name and additional information

    :param metaxml: file mit information for netcdf files (metadata)
    :param name: file name information
    :return: List with metadata information: metaNetcdfVar
    """
    if os.path.isfile(metaxml):
        try:
            metaparse = xml.dom.minidom.parse(metaxml)
        except:
            msg = "Error 303: using option file: " + metaxml
            raise CWATMError(msg)
    else:
        msg = "Cannot find option file: " + metaxml +"\n"
        path, name = os.path.split(metaxml)
        #metaxml = os.path.join(os.getcwd(), name)
        # using program name
        metaxml = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'cwatm', name)
        if os.path.isfile(metaxml):
            msg += "Using file: " + metaxml + " instead."
            print(CWATMWarning(msg))
        else:
            msg = "Error 304: Cannot find alternative option file: " + metaxml
            raise CWATMFileError(metaxml, msg, sname = name)

        try:
            metaparse = xml.dom.minidom.parse(metaxml)
        except:
            msg = "Error 305: Error using alternative option file: " + metaxml
            raise CWATMError(msg)

    # running through all output variable
    # if an output variable is not defined here the standard metadata is used
    # unit = "undefined", standard name = long name = variable name
    meta = metaparse.getElementsByTagName("CWATM")[0]

    for metavar in meta.getElementsByTagName("metanetcdf"):
        d = {}
        for key in list(metavar.attributes.keys()):
            if key != 'varname':
                d[key] = metavar.attributes[key].value
        key = metavar.attributes['varname'].value
        metaNetcdfVar[key] = d



