#!/bin/bash
#SBATCH --job-name=orac_lut_build        # Descriptive name for the job
#SBATCH --output=orac_lut_%j.log         # Output log file (%j = job ID)
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=1                # Number of CPU cores
#SBATCH --mem=4G                          # Memory allocation
#SBATCH --time=01:00:00                   # Maximum runtime
#SBATCH --partition=shared               # SLURM partition - change if needed

# Load environment variables and paths
source ~/.bash_profile


# monochromatic calculation
# aerosol section


export CREATE_ORAC_AEROSOL_LUT_DRIVER="input_files/driver/meteosat-11_seviri_aerosol.driver"
echo "making ..." $CREATE_ORAC_LUT_DRIVER
idl -e "create_orac_aerosol_lut_wrapper,srf_quad=1,mmfile='aerosol_a77.mm',lutfile='aerosol.lut',tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version=21"

export CREATE_ORAC_AEROSOL_LUT_DRIVER="input_files/driver/meteosat-11_seviri_aerosol.driver"
echo "making ..." $CREATE_ORAC_LUT_DRIVER
idl -e "create_orac_aerosol_lut_wrapper,srf_quad=1,mmfile='aerosol_a78.mm',lutfile='aerosol.lut',tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version=21"

export CREATE_ORAC_AEROSOL_LUT_DRIVER="input_files/driver/meteosat-11_seviri_aerosol.driver"
echo "making ..." $CREATE_ORAC_LUT_DRIVER
idl -e "create_orac_aerosol_lut_wrapper,srf_quad=1,mmfile='aerosol_a79.mm',lutfile='aerosol.lut',tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version=21"

