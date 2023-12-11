import sys
import numpy as np
import pandas as pd
import math
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy as apy
import astroplan as apl
from astropy.time import Time
from astropy.time import TimeDelta

# if you haven't already: pip install requests
import requests

class Request():

    def __init__(self, simbad_name):
        self.simbad_name = simbad_name
        self.goodFormat = False

        # Alternate names
        self.tic = None
        self.gaia_name = None

        # Coordinate info
        self.RA = None
        self.Dec = None
        self.pmRA = None
        self.pmDec = None
        self.epoch = None

        # Host star info
        self.Jmag = None
        self.Teff = None

        #Observation info
        self.nominal_ExpTime = None
        self.max_ExpTime = None
        self.simulcal = None
        self.eventWindow = None

        #Cadence info
        self.n_observations_per_visit = None
        self.n_visits_per_night = None
        self.n_unique_nights_per_semester = None
        self.intra_night_cadence = None
        self.inter_night_cadence = None

        #stats
        self.total_observations_requested = None

    def get_single_shot(self):
        self.get_star_info()
        #Cadence info
        self.n_observations_per_visit = 1
        self.n_visits_per_night = 1
        self.n_unique_nights_per_semester = 1
        self.intra_night_cadence = 0.0
        self.inter_night_cadence = 0

    def get_RM(self):
        self.get_star_info()
        #Cadence info
        self.n_observations_per_visit = math.ceil((self.eventWindow*3600)/self.nominal_ExpTime)
        self.n_visits_per_night = 1
        self.n_unique_nights_per_semester = 1
        self.intra_night_cadence = 0.0
        self.inter_night_cadence = 0

    # Function to retrieve basic coordinate info from SIMBAD
    def get_star_info(self):
        base_url = "http://simbad.u-strasbg.fr/simbad/sim-id"
        params = {
            'Ident': self.simbad_name,
            'output.format': 'ASCII'
        }

        response = requests.get(base_url, params=params)

        if response.status_code == 200:
            lines = response.text.split('\n')

            ra, dec, pm_ra, pm_dec, pm_epoch, ticname, gaianame, gaianame2, gaianame3 = None, None, None, None, None, None, None,  None, None

            identifiers_section = False
            bibcodes_section = False
            all_names = []

            j_mag = None
            for line in lines:
                if line.startswith('Bib'):
                    bibcodes_section = True
                if line.startswith('Identifiers ('):
                    identifiers_section = True
                if line.startswith('Coordinates(ICRS'):
                    coords = line.split(': ')[1].split()
                    ra, dec = ' '.join(coords[:3]), ' '.join(coords[3:])
                    pm_epoch = line.split('=')[1].strip()[:5]

                elif line.startswith('Flux J'):
                    j_mag = float(line.split(' ')[3].strip())
                elif line.startswith('Proper motions'):
                    pm_values = line.split(':')[1].split()[:2]
                    pm_ra, pm_dec = float(pm_values[0]), float(pm_values[1])
                elif identifiers_section == True and bibcodes_section == False:
                    identifiers = line.split()
                    for i in identifiers:
                        all_names.append(i)

            for n in range(len(all_names)):
                if all_names[n] == 'TIC':
                    ticname = 'TIC' + str(all_names[n+1])
                if all_names[n] == 'Gaia' and all_names[n+1] == 'DR1':
                    gaianame = 'Gaia_DR1_' + str(all_names[n+2])
                if all_names[n] == 'Gaia' and all_names[n+1] == 'DR2':
                    gaianame2 = 'Gaia_DR2_' + str(all_names[n+2])
                if all_names[n] == 'Gaia' and all_names[n+1] == 'DR3':
                    gaianame3 = 'Gaia_DR3_' + str(all_names[n+2])

            finalgaianame = None
            if gaianame != None:
                finalgaianame = gaianame
            else:
                if gaianame2 != None:
                    finalgaianame = gaianame2
                else:
                    if gaianame3 != None:
                        finalgaianame = gaianame3

            if ra is None or dec is None:
                print("Star not found. Error 1.")
                return None
            else:
                rasplit = ra.split(" ")
                ra = rasplit[0] + ":" + rasplit[1] + ":" + rasplit[2][:5]
                decsplit = dec.split(" ")
                dec = decsplit[0] + ":" + decsplit[1] + ":" + decsplit[2][:5]

                self.tic = ticname
                self.gaia_name = finalgaianame
                self.RA = ra
                self.Dec = dec
                self.pmRA = pm_ra
                self.pmDec = pm_dec
                self.epoch = pm_epoch
                self.Jmag = j_mag
                # Convert SIMBAD default RA and Dec format to decimal format
                c = SkyCoord(ra=self.RA, dec=self.Dec, unit=(u.hourangle, u.deg))
                ra_angle = c.ra.deg
                dec_angle = c.dec.deg
                self.RA = round(ra_angle,2)
                self.Dec = round(dec_angle,2)
        else:
            print("Star not found. Error 2.")


    def computeTimeRequest(self):
        # Do not change these values.
        readout_time = 45 #seconds
        avg_slew = 240 #seconds

        self.total_observations_requested = self.n_observations_per_visit*self.n_visits_per_night*self.n_unique_nights_per_semester
        self.overhead_per_unique_night = (self.n_observations_per_visit-1)*readout_time + (self.n_visits_per_night*avg_slew)
        self.time_per_unique_night = self.nominal_ExpTime + self.overhead_per_unique_night
        self.total_time_for_target_seconds = self.time_per_unique_night*self.total_observations_requested
        self.total_time_for_target_hours = round(self.total_time_for_target_seconds/3600,2)


    def isObservable(self):

        min_az = 5.3
        max_az = 146.2
        min_alt = 33.3
        else_min_alt = 25
        max_alt = 85

        twilight_frame = pd.read_csv("twilight_times_2024A.csv")
        coords = apy.coordinates.SkyCoord(self.RA * u.hourangle, self.Dec * u.deg, frame='icrs')
        target = apl.FixedTarget(name=self.simbad_name, coord=coords)

        semesterly_observability_matrix = []
        for i in range(len(twilight_frame)):

            start_time = twilight_frame['12_evening'][i]
            start_time_jd = Time(start_time,format='jd')
            end_time = twilight_frame['12_morning'][i]
            end_time_jd = Time(end_time,format='jd')

            # test in steps of 5 minutes
            step = TimeDelta(300,format='sec')
            nightly_steps = np.arange(start_time_jd.jd, end_time_jd.jd, step.jd)
            t = Time(nightly_steps,format='jd')

            keck = apl.Observer.at_site('W. M. Keck Observatory')
            AZ = keck.altaz(t, target, grid_times_targets=True)

            for j in range(len(AZ)):
                alt=AZ[j].alt.deg
                az=AZ[j].az.deg
                deck = np.where((az >= min_az) & (az <= max_az))
                deck_height = np.where((alt <= max_alt) & (alt >= min_alt))
                first = np.intersect1d(deck,deck_height)

                not_deck_1 = np.where((az < min_az))
                not_deck_2 = np.where((az > max_az))
                not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt))
                second = np.intersect1d(not_deck_1,not_deck_height)
                third = np.intersect1d(not_deck_2,not_deck_height)

                good = np.concatenate((first,second,third))

                nightly_observability_matrix = np.zeros(len(t),dtype=int)
                nightly_observability_matrix[good] = 1

                # test if up for more than 20 minutes (~1/8 of a quarter night)
                if np.sum(nightly_observability_matrix) > 4:
                    semesterly_observability_matrix.append(1)
                else:
                    semesterly_observability_matrix.append(0)

        self.observability = semesterly_observability_matrix
        if np.sum(self.observability) >= 1:
            self.first_observable = semesterly_observability_matrix.index(1)
            self.last_observable = len(semesterly_observability_matrix) - 1 - semesterly_observability_matrix[::-1].index(1)
            self.days_observable = self.last_observable - self.first_observable
            print("This target is observable for a total of [" + str(self.days_observable) + "] days this semester.")
            print("This target rises on day [" + str(self.first_observable) + "] and sets on day [" + str(self.last_observable) + "] of the semester. \n")
            if self.days_observable <= 30:
                print("Further warning: this target is hardly accessible from Hawaii (accessible for less than 30 days this semester).")
        else:
            print("WARNING: this target is not accessible from Hawaii at all.")
            self.first_observable = 0
            self.last_observable = 0
            self.days_observable = 0

    def determineFeasibility(self):

        historical_unique_nights_allocated = 60
        historical_quarters_allocated = 120

        if self.n_unique_nights_per_semester > historical_unique_nights_allocated:
            print("Your program may not be feasibile. Historically, HIRES/KPF is allocated " + str(historical_unique_nights_allocated) + " unique nights in a semester. You have asked for a greater number of unique night visits for your request. Consider ammending your request, or accept that it may not be able to be fully completed. \n")
        elif self.n_unique_nights_per_semester*self.inter_night_cadence > 180:
            print("Your program may not be feasible. With your desired inter-night cadence and desired unique night visits, we would require more than 180 days in the semester to achieve your request. Consider ammending your request, or accept that it may not be able to be fully completed. \n")
        elif self.n_unique_nights_per_semester*self.inter_night_cadence > self.days_observable:
            print("Your program may not be feasible. With your desired inter-night cadence and desired unique night visits, combined with this target's accessibility this semester, there is likely not enough unique nights in the semester to complete this request. Consider ammending your request, or accept that it may not be able to be fully completed. \n")
        elif self.days_observable < historical_unique_nights_allocated:
            print("Your program may not be feasible. Your target is observable for fewer unique nights than you have requested for observations. Consider ammending your request, or accept that it may not be able to be fully completed. \n")
        else:
            print("To first order checks, this target's request is feasible. This is not a guarentee that it will actually be completed to 100% of the request.")

    def runChecks(self):

        # First check for Nones
        print("Checking for parameters still set to None.")
        allattributes = {
                         "TIC":self.tic,
                         "Gaia name":self.gaia_name,
                         "RA":self.RA,
                         "Dec":self.Dec,
                         "pmRA":self.pmRA,
                         "pmDec":self.pmDec,
                         "Epoch":self.epoch,
                         "Jmag":self.Jmag,
                         "Teff":self.Teff,
                         "Nominal ExpTime":self.nominal_ExpTime,
                         "Max ExpTime":self.max_ExpTime,
                         "SimulCal":self.simulcal,
                         "N_observations_per_visit":self.n_observations_per_visit,
                         "N_visits_per_night":self.n_visits_per_night,
                         "N_unique_nights_per_semester":self.n_unique_nights_per_semester,
                         "Intra-night cadence":self.intra_night_cadence,
                         "Inter-night cadence":self.inter_night_cadence
                        }

        keys = allattributes.keys()
        self.notNones = True
        for k in keys:
            if allattributes[k] == None:
                print(str(k) + " is None.")
                self.notNones = False

        print()
        print("Checking for correct data types.")
        self.goodTypes = True

        if isinstance(self.n_observations_per_visit, int) == False:
            print("N_observations_per_visit must be an integer. Fix before continuing")
            self.goodTypes = False
        if isinstance(self.n_visits_per_night, int) == False:
            print("N_visits_per_night must be an integer. Fix before continuing")
            self.goodTypes = False
        if isinstance(self.n_unique_nights_per_semester, int) == False:
            print("N_unique_nights_per_semester must be an integer. Fix before continuing")
            self.goodTypes = False
        if isinstance(self.intra_night_cadence, float) == False:
            print("Intra_night_cadence must be a float. Fix before continuing")
            self.goodTypes = False
        if isinstance(self.inter_night_cadence, int) == False:
            print("Inter_night_cadence must be an integer. Fix before continuing")
            self.goodTypes = False

        print()
        print("Checking for formatting.")
        self.goodFormats = True

        if self.n_visits_per_night == 1 and self.intra_night_cadence != 0:
            print("Only one visit per night requested but intra-night cadence is greater than 0. Setting intra-night cadence to 0 (default).")
            self.intra_night_cadence = 0
        if self.max_ExpTime < 0:
            print("Max exposure time must be greater than or equal to nominal exposure time. Resetting to be equal. \n")
            self.max_ExpTime = self.nominal_ExpTime
        if self.nominal_ExpTime <= 0:
            print("No negative nominal exposure times. Please retry.\n")
            self.goodFormats = False
        if self.n_observations_per_visit < 1:
            print("Minimum observations per visit is 1. Please retry.\n")
            self.goodFormats = False
        if self.intra_night_cadence < 0.5 and self.intra_night_cadence != 0:
            print("Minimum intra-night cadence is 0.5 hours. Please retry.\n")
            self.goodFormats = False
        if self.n_unique_nights_per_semester < 1:
            print("Minimum unique nights per semester is 1. Please retry.\n")
            self.goodFormats = False
        print()

        if self.goodTypes == False:
            print("At least one of your data types is incorrect. Fix before continuing.")
        if self.goodFormats == False:
            print("At least one of your data formats is incorrect. Fix before continuing.")
        if self.notNones == False:
            print("At least one of your parameters is None. Fix before continuing.")
        if self.goodTypes and self.goodFormats and self.notNones:
            print("All looks good for this request! \n")
            self.canContinue = True
        else:
            self.canContinue = False

class Program():

    def __init__(self, requests):
        self.requests = requests

        # admin
        self.code = None
        self.semester = None
        self.savefile = None

        # stats
        self.totalProgramObservations = None
        self.total_time_for_program_seconds = None
        self.total_time_for_program_hours = None
        self.total_time_for_program_nights = None

    def totalProgramTime(self):

        # Do not change these values.
        readout_time = 45 #seconds
        avg_slew = 240 #seconds

        self.totalProgramObservations = 0
        self.total_time_for_program_seconds = 0

        for i in range(len(self.requests)):
            self.totalProgramObservations += self.requests[i].total_observations_requested
            self.total_time_for_program_seconds += self.requests[i].total_time_for_target_seconds
        self.total_time_for_program_hours = round(self.total_time_for_program_seconds/3600,2)
        self.total_time_for_program_nights = round((self.total_time_for_program_seconds/3600)/10,2)
        print("This program requires an allocation of [" + str(self.total_time_for_program_nights) + "] nights to be feasible.")

    def checkAllRequests(self):

        allGood = True
        for i in range(len(self.requests)):
            if self.requests[i].goodTypes == False or self.requests[i].goodFormats == False or self.requests[i].notNones == False:
                allGood = False
        self.allGood = allGood

    def writeFile(self):

        if self.allGood:

            with open(self.savefile,'w') as fileout:
                fileout.write("Name,Gaia_ID,TIC_ID,Program_Code,RA,Dec,pmRA,pmDec,Epoch,Jmag,Teff,Nominal_ExpTime,Max_ExpTime,Simulcal,Total_Observations_Requested,N_Observations_per_Visit,N_Visits_per_Night,Intra_Night_Cadence,N_Unique_Nights,Inter_Night_Cadence,Total_Time_for_Target_Seconds,Total_Time_for_Target_Hours" + "\n")

                for r in range(len(self.requests)):

                    line = str(self.requests[r].simbad_name) + "," + str(self.requests[r].gaia_name) + "," + str(self.requests[r].tic) + "," + str(self.code) + "," + \
                    str(self.requests[r].RA) + "," + str(self.requests[r].Dec) + "," + \
                    str(self.requests[r].pmRA) + "," + str(self.requests[r].pmDec) + "," + \
                    str(self.requests[r].epoch) + "," + \
                    str(self.requests[r].Jmag) + "," + str(self.requests[r].Teff) + "," + \
                    str(self.requests[r].nominal_ExpTime) + "," + str(self.requests[r].max_ExpTime) + "," + \
                    str(self.requests[r].simulcal) + "," + \
                    str(self.requests[r].total_observations_requested) + "," + str(self.requests[r].n_observations_per_visit) + "," + \
                    str(self.requests[r].n_visits_per_night) + "," + str(self.requests[r].intra_night_cadence) + "," + \
                    str(self.requests[r].n_unique_nights_per_semester) + "," + str(self.requests[r].inter_night_cadence) + "," + \
                    str(self.requests[r].total_time_for_target_seconds) + "," + str(self.requests[r].total_time_for_target_hours)

                    print(line)
                    print()
                    fileout.write(line + "\n")

        else:
            print("Not all targets are formatted correctly. Will not write file until all checks pass.")
