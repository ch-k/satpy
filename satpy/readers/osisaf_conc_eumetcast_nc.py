#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""An OSISAF SST reader for the netCDF GHRSST format
"""

import logging
from datetime import datetime
import numpy as np
import xarray as xr

from satpy.dataset import Dataset
from satpy.readers.netcdf_utils import NetCDF4FileHandler

logger = logging.getLogger(__name__)

PLATFORM_NAME = {'NPP': 'Suomi-NPP', }
SENSOR_NAME = {'VIIRS': 'viirs',
               'AVHRR': 'avhrr/3'}


class OSISAF_CONC_EUMETCAST_NC(NetCDF4FileHandler):
    """OSISAF ice concentration EUMETCast netCDF reader."""

    def _parse_datetime(self, datestr):
        return datetime.strptime(datestr, '%Y-%m-%d %H:%M:%S')

    def get_area_def(self, area_id):
        """Override abstract baseclass method"""
        from pyresample.geometry import AreaDefinition
        from pyresample import create_area_def
    #     from pyresample.utils import proj4_str_to_dict
        proj4str = self['Polar_Stereographic_Grid/attr/proj4_string']
        x_size = self['/dimension/xc']
        y_size = self['/dimension/yc']
        p_lowerleft_lat = self['lat'].values[y_size-1,0]
        p_lowerleft_lon = self['lon'].values[y_size-1,0]
        p_upperright_lat = self['lat'].values[0,x_size-1]
        p_upperright_lon = self['lon'].values[0,x_size-1]
        area_extent = [p_lowerleft_lon, p_lowerleft_lat, p_upperright_lon, p_upperright_lat]
        area_def = create_area_def(area_id='osisaf_conc_polar_stereographic', 
            description='osisaf_conc_polar_stereographic',
            proj_id='osisaf_conc_polar_stereographic',
            projection=proj4str, width=x_size, height=y_size, area_extent=area_extent, units='deg')
        return area_def
        # raise NotImplementedError

    def get_dataset(self, dataset_id, ds_info):
        """Load a dataset."""

        logger.debug('Reading {} from {}'.format(dataset_id.name, self.filename))
        var_path = ds_info.get('file_key', '{}'.format(dataset_id.name))
        dtype = ds_info.get('dtype', np.float32)
        if var_path + '/shape' not in self:
            # loading a scalar value
            shape = 1
        else:
            shape = self[var_path + '/shape']
            if shape[0] == 1:
                # Remove the time dimenstion from dataset
                shape = shape[1], shape[2]

        file_units = ds_info.get('file_units')
        if file_units is None:
            try:
                file_units = self[var_path + '/attr/units']
                # they were almost completely CF compliant...
                if file_units == "none":
                    file_units = "1"
            except KeyError:
                # no file units specified
                file_units = None

        data = self[var_path][0]
        valid_min = self[var_path + '/attr/valid_min']
        valid_max = self[var_path + '/attr/valid_max']
        if valid_min is not None and valid_max is not None:
            data = data.where(data >= valid_min, np.nan)
            data = data.where(data <= valid_max, np.nan)
        try:
            scale_factor = self[var_path + '/attr/scale_factor']
            scale_offset = self[var_path + '/attr/add_offset']
        except KeyError:
            scale_factor = scale_offset = None

        data = (data * scale_factor + scale_offset)
        
        # Set proper dimension names
        data = data.rename({'xc': 'x', 'yc': 'y'})

        ds_info.update({
            "units": ds_info.get("units", file_units),
            "platform_name": PLATFORM_NAME.get(self['/attr/platform_name'], self['/attr/platform_name']),
            "sensor": SENSOR_NAME.get(self['/attr/instrument_type'], self['/attr/instrument_type']),
        })
        ds_info.update(dataset_id.to_dict())
        data.attrs.update(ds_info)
        return data

    # def get_lonlats(self, navid, nav_info, lon_out=None, lat_out=None):
    #     """Load an area.
    #     """
    #     lon_key = 'lon'
    #     valid_min = self[lon_key + '/attr/valid_min']
    #     valid_max = self[lon_key + '/attr/valid_max']

    #     lon_out.data[:] = self[lon_key][::-1]
    #     lon_out.mask[:] = (lon_out < valid_min) | (lon_out > valid_max)

    #     lat_key = 'lat'
    #     valid_min = self[lat_key + '/attr/valid_min']
    #     valid_max = self[lat_key + '/attr/valid_max']
    #     lat_out.data[:] = self[lat_key][::-1]
    #     lat_out.mask[:] = (lat_out < valid_min) | (lat_out > valid_max)

    #     return {}

    @property
    def start_time(self):
        return self.filename_info['start_time']
        # return self._parse_datetime(self['/attr/start_date'])

    @property
    def end_time(self):
        return self._parse_datetime(self['/attr/stop_date'])
