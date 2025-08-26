from dataclasses import dataclass
from os.path import isfile
import re
from ..utils.external import *
# The `make_calc_difference_sf` and `make_difference_map` functions below is related to the `sfallFoFo.sh` of the original bash script

@dataclass
class CCP4_DED_Utils:
    col_names: dict[str, str | None]
    high_resolution: float
    low_resolution: float
    symmetry: str
    grid_ccp4: str
    default_sigma_level: float
    atom_pattern: list[str]
    should_use_phenix_for_dFc: bool = False

    def make_calc_difference_sf(self, model_in: str, data_out: str, dark_data: str, light_data: str | None, dark_model_in: str | None = None, keep_files=False):
        """@param model_in: input structure model file ending in `.pdb`"""
        # structure factor calculation
        ctrl=f'''TITL dark_lit_model
                MODE SFCALC XYZIN HKLIN
                LABI FP={self.col_names['Fo_dark']} SIGFP={self.col_names['Fo_dark_sigma']}
                LABO FC=FCalclight PHIC=PHICalclight
                RESO {self.low_resolution} {self.high_resolution}
                BINS 60
                END
            '''
        log = sfall(
            xyzin=model_in,
            hklin=dark_data,
            hklout=(data_of_model := model_in.replace('.pdb', '.mtz')),
            script=ctrl
        )
        log += u'\u2500' * 30 + '\n'

        # calculate difference structure factors between model and dark
        if self.should_use_phenix_for_dFc:
            if not dark_model_in:
                raise RuntimeError('Dark model should be given as the phase source of `phenix.fobs_minus_fobs_map`')
            log += fobs_minus_fobs_map(
                job_title="calculate_dFc_between_model_and_dark",
                f_obs_1_file=f'{data_of_model}',
                f_obs_2_file=f'{dark_data}',
                output_file=f'{data_out}',
                f_obs_1_label="FCalclight",
                f_obs_2_label="FCalcdark",
                phase_source=f'{dark_model_in}',
                high_res=self.high_resolution,
                low_res=self.low_resolution,
                sigma_cutoff=self.default_sigma_level,
            )
        else:
            ctrl=f'''read {dark_data}
                    read {data_of_model}
                    {f"read {light_data}" if light_data and light_data.endswith('.mtz') else ""}
                    calc F col dFc = col FCalclight col {self.col_names['Fc_dark']} -
                    calc P col corr = 0 DTR
                    select col dFc < 0
                    calc P col corr = 180 DTR
                    calc F col dFc = col dFc -1 *
                    select all
                    calc P col dPHIc = col {self.col_names['Fc_dark_phase']} col corr +
                    absent col dFc dPHIc if  col FoFo = absent
                    write {data_out} col dFc dPHIc
                    END
                '''
            log += sftools(script=ctrl)

        # remove all used files unless we need them (e.g. regenerate after each round)
        used_files = [
            data_of_model, # structure factor of the generated model
        ]
        if not keep_files:
            for used_file in used_files:
                subprocess.run(['rm', used_file])

        # Checkpoint: check if the last output file exists 
        if not isfile(data_out):
            with open(f"log_dFoCC.make_calc_difference_sf.txt", "w") as f:
                f.write(log)
            raise Exception(log)
        
        return log

    def make_difference_map(self, data_in: str, map_out: str, is_dFo=False):
        # export map with the same characteristics as the observed difference structure factors
        is_from_phenix = self.should_use_phenix_for_dFc or is_dFo # FIXME: rename is_dFo to is_from_phenix
        ctrl=f'''RESOLUTION {self.low_resolution} {self.high_resolution}
                LABIN F1={"FoFo" if is_from_phenix else "dFc"} PHI={"PHFc" if is_from_phenix else "dPHIc"}
                {'' if self.grid_ccp4 == 'SAMPLE 2' else f"GRID {self.grid_ccp4}"}
                END
            '''
        log = fft(
            hklin=data_in,
            mapout=map_out,
            script=ctrl
        )

        # Checkpoint: check if the last output file exists 
        if not isfile(map_out):
            with open(f"log_dFoCC.make_difference_map.txt", "w") as f:
                f.write(log)
            raise Exception(log)
        
        return log

    # The `calculate_sigma_from_mean_density` and `generate_sigmadiff_maps` functions below are related to the `correlate_all.sh` of the original bash script
    def calculate_sigma_from_mean_density(self, map_in: str):
        # get the value of sigma for the entire map
        ctrl='''GO
                END
            '''
        log = mapdump(
            # mapin=f'dFc.{prefix}_above0.map',
            mapin=map_in,
            script=ctrl
        )

        for line in log.splitlines():
            match = re.search(r'Rms deviation from mean density...................(....\d......)', line)
            if match:
                rms = float(match.group(1))

        # Checkpoint: check if the last output value is defined
        try:
            rms
        except:
            with open("log_dFoCC.calculate_sigma_from_mean_density.txt", "w") as f:
                f.write(log)
            raise Exception(log)

        return log, rms

    def generate_sigmadiff_maps(self, map_in: str, map_out: str, rms: float, sigma_level: float, keep_files=False):
        """@param map_in: input map file ending in `.map`"""
        sigma_level = sigma_level if sigma_level else self.default_sigma_level
        sigma = rms * sigma_level
        
        # make a positive sigma mask
        ctrl=f'''SYMMETRY {self.symmetry}
                MASK -
                CUT {sigma}
                END
            '''
        log = mapmask(
            mapin=map_in,
            mskout=(pos_sigma_mask := map_in.replace('.map', f'_above{self.default_sigma_level}.msk')),
            script=ctrl
        )
        
        # make a negative sigma mask
        ctrl=f'''SYMMETRY {self.symmetry}
                MASK -
                CUT -{sigma}
                END
            '''
        log += mapmask(
            mapin=map_in,
            mskout=(neg_sigma_mask := map_in.replace('.map', f'_below{self.default_sigma_level}.msk')),
            script=ctrl
        )
        
        # make a positive map (remove signals below sigma)
        ctrl='''MAP INCLUDE
                END
            '''
        log += overlapmap(
            mapin1=pos_sigma_mask,
            mapin2=map_in,
            mapout=(pos_sigma_map := map_in.replace('.map', f'_all_above{self.default_sigma_level}.map')),
            script=ctrl
        )
        
        # make a negative map (remove signals above negative sigma)
        ctrl='''MAP EXCLUDE
                END
            '''
        log += overlapmap(
            mapin1=neg_sigma_mask,
            mapin2=map_in,
            mapout=(neg_sigma_map := map_in.replace('.map', f'_all_below{self.default_sigma_level}.map')),
            script=ctrl
        )
        
        # combine positive and negative maps to get +/- sigma cleared map
        ctrl='''MAP ADD
                END
            '''
        log += overlapmap(
            mapin1=pos_sigma_map,
            mapin2=neg_sigma_map,
            mapout=map_out,
            script=ctrl
        )
        
        # remove all used files
        used_files = [
            # map_in, # useful in certain cases
            pos_sigma_mask,
            neg_sigma_mask,
            pos_sigma_map,
            neg_sigma_map,
            # map_out, # still in use
        ]
        # keep the calculated difference map if we need it
        if not keep_files:
            for used_file in used_files:
                subprocess.run(['rm', used_file])

        # Checkpoint: check if the last output file exists
        if not isfile(map_out):
            with open(f"log_dFoCC.generate_sigmadiff_maps.txt", "w") as f:
                f.write(log)
            raise Exception(log)

        return log


    # # The `generate_sigma_cutoff_map` replaces the `calculate_sigma_from_mean_density` and `generate_sigmadiff_maps` functions
    # # 
    # # It uses a library called `gemmi` instead of CCP4
    # def generate_sigma_cutoff_map(self, map_in: str, map_out: str, sigma_level: float):
    #     """
    #     Cut off the noise density between +/- certain sigma level
    #     @param sigma_level: default is 3.0 sigma level 
    #     """
    #     sigma_level = sigma_level if sigma_level else self.default_sigma_level
    #
    #     ccp4_map = gemmi.read_ccp4_map(map_in, setup=True)
    #     sigma_value = ccp4_map.grid.array.std()
    #     cutoff = sigma_value * sigma_level
    #     ccp4_map.grid.array[(ccp4_map.grid.array >= -cutoff) & (ccp4_map.grid.array <= +cutoff)] = 0.
    #     ccp4_map.update_ccp4_header(update_stats=True)
    #     ccp4_map.write_ccp4_map(map_out)
    #
    #     # # remove all used files
    #     # used_files = [
    #     #     # f'dFc.{prefix}.map' # useful in certain cases
    #     #     # f'{workdir}dFc.{prefix}_sigma{sigma_level}.map', # still in use
    #     # ]
    #     # # keep the calculated difference map if we need it
    #     # if not keep_files:
    #     #     for used_file in used_files:
    #     #         subprocess.run(['rm', used_file])

    # The `extract_region` function below is related to `sfall6ns.sh` of the original bash script
    # It grabs only the atoms in the moiety, and make a new pdb for masking
    # 
    # If anything about the moiety is changed, this part must be modified
    def extract_region(self, xyzin_ref: str, xyzin: str, xyzout: str):
        subprocess.run(f"grep -v 'END' {xyzin_ref} > {xyzout}", shell=True)
        for atom_pattern in self.atom_pattern:
            subprocess.run([f"grep -i '{atom_pattern}' {xyzin} >> {xyzout}"], shell=True)
        subprocess.run([f'echo "END" >> {xyzout}'], shell=True)

        # subprocess.run(['cp', xyzin_ref, xyzout])
        # subprocess.run([f"grep -i 'N   LYS A 216' {xyzin} >> {xyzout}"], shell=True)
        # # ...

        # Checkpoint: check if the last output file exists
        if not isfile(xyzout):
            log = f"{xyzin_ref = }\n{xyzin = }\n{xyzout = }"
            with open(f"log_dFoCC.extract_region.txt", "w") as f:
                f.write(log)
            raise Exception(log)
