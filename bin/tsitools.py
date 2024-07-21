"""
Author: Tayeb Kakeshpour
Description: Tools for analyzing TSI 3330 output files
"""
import pandas as pd
import math,os
import numpy as np
import datetime

class ParticleData:
    def __init__(self, directory, file_suffix, num_time_points=None):
        self.directory = directory
        self.file_suffix = file_suffix
        self.num_time_points = num_time_points
        self.cut_points = {}
        self.mean_diameters = {}
        self.particle_volumes = {}
        self.times = []
        self.particle_counts = {}

        file_path = self._find_file()
        if file_path:
            self.file_path = file_path
            self._parse_file()
            self._calculate_mean_diameters()
            self._calculate_particle_volumes()
        else:
            print("No file found with the specified suffix in the directory.")
        
        self._test_start_date = self._extract_test_start_date()
        self._test_start_time = self._extract_test_start_time()

    def _find_file(self):
        for file in os.listdir(self.directory):
            if file.endswith(self.file_suffix + '.csv'):
                return os.path.join(self.directory, file)
        return None

    def _parse_file(self):
        # Reading the file into a list of lines
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Extracting the cut points
        for line in lines:
            if line.startswith('Bin '):
                parts = line.split(',')
                bin_number = parts[0].split(' ')[0] + ' ' + parts[0].split(' ')[1]  # Extract 'Bin X'
                cut_point = float(parts[1].split()[0])
                self.cut_points[bin_number] = cut_point

        # Finding the start of the data rows
        data_start_line = next(i for i, line in enumerate(lines) if line.startswith('Elapsed Time [s]'))

        # Reading the data rows into a DataFrame
        data = pd.read_csv(self.file_path, skiprows=data_start_line)
        
        if self.num_time_points is not None:
            data = data.head(self.num_time_points)
        
        # Extracting the time and particle count data
        self.times = data['Elapsed Time [s]'].tolist()
        for column in data.columns:
            if column.startswith('Bin '):
                self.particle_counts[column] = data[column].tolist()

    def _extract_test_start_date(self):
        # Reading the file into a list of lines
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Extracting the test start date
        for line in lines:
            if line.startswith('Test Start Date'):
                date_str = line.split(',')[1].strip()
                date_obj = datetime.datetime.strptime(date_str, '%Y/%m/%d')
                formatted_date = date_obj.strftime('%m/%d/%y')  # American format MM/DD/YY
                return formatted_date
        return None

    def _extract_test_start_time(self):
        # Reading the file into a list of lines
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Extracting the test start time
        for line in lines:
            if line.startswith('Test Start Time'):
                time_str = line.split(',')[1].strip()
                time_obj = datetime.datetime.strptime(time_str, '%H:%M:%S')
                formatted_time = time_obj.strftime('%H:%M')  # Format HH:MM
                return formatted_time
        return None
    
    def _calculate_mean_diameters(self):
        for i in range(1, 17):
            lower_cut_point = self.cut_points[f"Bin {i}"]
            upper_cut_point = self.cut_points[f"Bin {i+1}"]
            self.mean_diameters[f"Bin {i}"] = (lower_cut_point + upper_cut_point) / 2

        # Last bin (Bin 17) uses 10.3 micron
        self.mean_diameters["Bin 17"] = 10.3

    def _calculate_dlogD(self):
        self.dlogD = {}
        for i in range(1, 17):
            lower_cut_point = self.cut_points[f"Bin {i}"]
            upper_cut_point = self.cut_points[f"Bin {i+1}"]
            self.dlogD[f"Bin {i}"] = np.log10(upper_cut_point) - np.log10(lower_cut_point)
        self.dlogD["Bin 17"] = 0.095

    def _calculate_particle_volumes(self):
        for bin_name, diameter in self.mean_diameters.items():
            radius = diameter / 2
            volume = (4/3) * math.pi * (radius ** 3)
            self.particle_volumes[bin_name] = volume

    def get_dlogD_dict(self):
        if not hasattr(self, 'dlogD'):
            self._calculate_dlogD()
        return self.dlogD

    def get_dlogD_np(self):
        dlogD_dict = self.get_dlogD_dict()
        return np.array(list(dlogD_dict.values()))

    def get_cut_points(self):
        return self.cut_points

    def get_mean_diameters(self):
        return self.mean_diameters

    def get_particle_volumes(self):
        return self.particle_volumes

    def get_times(self):
        return self.times

    def get_particle_counts(self, bin_number):
        return self.particle_counts.get(f'Bin {bin_number}', None)

    def get_particle_counts_at_time(self, time):
        if time not in self.times:
            return None  # Time not found
        index = self.times.index(time)
        counts_at_time = {bin_number: counts[index] for bin_number, counts in self.particle_counts.items()}
        return counts_at_time

    def get_average_particle_counts(self):
        averages = {}
        for bin_name, counts in self.particle_counts.items():
            averages[bin_name] = np.mean(counts)
        return averages

    def get_total_volume_per_bin(self, use_average_counts=True, specific_time=None):
        if use_average_counts:
            counts = self.get_average_particle_counts()
        else:
            if specific_time is None:
                raise ValueError("Specific time must be provided if not using average counts.")
            counts = self.get_particle_counts_at_time(specific_time)
        total_volumes = {bin_name: self.particle_volumes[bin_name] * count for bin_name, count in counts.items()}
        return total_volumes


    def get_total_particle_volume(self, use_average_counts=True, specific_time=None):
        if use_average_counts:
            counts = self.get_average_particle_counts()
        else:
            if specific_time is None:
                raise ValueError("Specific time must be provided if not using average counts.")
            counts = self.get_particle_counts_at_time(specific_time)

        total_volume = 0
        for bin_name, count in counts.items():
            volume_per_particle = self.particle_volumes.get(bin_name, 0)
            total_volume += volume_per_particle * count

        return total_volume

    def get_average_particle_counts_np(self):
        averages = self.get_average_particle_counts()
        return np.array(list(averages.values()))

    def get_total_volume_per_bin_np(self, use_average_counts=True, specific_time=None):
        total_volumes = self.get_total_volume_per_bin(use_average_counts, specific_time)
        return np.array(list(total_volumes.values()))
    



    
    def get_dN_dlogD(self, use_average_counts=True, specific_time=None):
        if use_average_counts:
            counts = self.get_average_particle_counts()
        else:
            if specific_time is None:
                raise ValueError("Specific time must be provided if not using average counts.")
            counts = self.get_particle_counts_at_time(specific_time)

        # Ensure dlogD is calculated
        if not hasattr(self, 'dlogD'):
            self._calculate_dlogD()

        dN_dlogD = {}
        for bin_name in counts.keys():
            count = counts[bin_name]
            dlogD_value = self.dlogD.get(bin_name, 1)  # Default to 1 to avoid division by zero
            dN_dlogD[bin_name] = count / dlogD_value if dlogD_value else 0

        return dN_dlogD



    def get_dN_dlogD_dict(self, use_average_counts=True, specific_time=None):
        dN_dlogD = self.get_dN_dlogD(use_average_counts, specific_time)
        return dN_dlogD

    def get_dN_dlogD_np(self, use_average_counts=True, specific_time=None):
        dN_dlogD_dict = self.get_dN_dlogD_dict(use_average_counts, specific_time)
        return np.array(list(dN_dlogD_dict.values()))

    def get_dV_dlogD(self, use_average_counts=True, specific_time=None):
        if use_average_counts:
            counts = self.get_average_particle_counts()
        else:
            if specific_time is None:
                raise ValueError("Specific time must be provided if not using average counts.")
            counts = self.get_particle_counts_at_time(specific_time)

        # Ensure dlogD is calculated
        if not hasattr(self, 'dlogD'):
            self._calculate_dlogD()

        dV_dlogD = {}
        for bin_name, count in counts.items():
            volume_per_particle = self.particle_volumes.get(bin_name, 0)
            volume_total = volume_per_particle * count
            dlogD_value = self.dlogD.get(bin_name, 1)  # Default to 1 to avoid division by zero
            dV_dlogD[bin_name] = volume_total / dlogD_value if dlogD_value else 0

        return dV_dlogD

    def get_dV_dlogD_dict(self, use_average_counts=True, specific_time=None):
        dV_dlogD = self.get_dV_dlogD(use_average_counts, specific_time)
        return dV_dlogD

    def get_dV_dlogD_np(self, use_average_counts=True, specific_time=None):
        dV_dlogD_dict = self.get_dV_dlogD_dict(use_average_counts, specific_time)
        return np.array(list(dV_dlogD_dict.values()))

    def get_mean_diameters_np(self):
        # Extracting mean diameters and converting them to a NumPy array
        mean_diameters_array = np.array(list(self.mean_diameters.values()))
        return mean_diameters_array

    def _calculate_bin_widths(self):
        self.bin_widths = {}
        cut_points_list = list(self.cut_points.values())

        for i in range(len(cut_points_list) - 1):
            bin_name = f"Bin {i+1}"
            self.bin_widths[bin_name] = cut_points_list[i+1] - cut_points_list[i]

        # Duplicate the width of the last bin
        last_bin_name = f"Bin {len(cut_points_list)}"
        self.bin_widths[last_bin_name] = self.bin_widths[f"Bin {len(cut_points_list)-1}"]

    def get_bin_widths(self):
        if not hasattr(self, 'bin_widths'):
            self._calculate_bin_widths()
        return self.bin_widths

    def get_bin_widths_np(self):
        bin_widths_dict = self.get_bin_widths()
        return np.array(list(bin_widths_dict.values()))
    
    def get_date(self):
        """Returns the test start date in MM/DD/YY format."""
        return self._test_start_date

    def get_time(self):
        """Returns the test start time in HH:MM format."""
        return self._test_start_time
    
    def get_particle_counts_at_time_np(self, time):
        """Get particle counts at a specific time as a NumPy array.

        Args:
            time (float): The time for which to get the particle counts.

        Returns:
            np.ndarray: The particle counts at the specified time, or None if the time is not found.
        """
        counts_at_time = self.get_particle_counts_at_time(time)
        if counts_at_time is None:
            return None  # Time not found, or some error occurred
        
        # Convert the dictionary values to a NumPy array
        counts_array = np.array(list(counts_at_time.values()))
        return counts_array

