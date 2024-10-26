using Distributed
using FITSIO
using Photometry
using DelimitedFiles
using WCS
using Reproject
using NaNStatistics

procs = addprocs(Sys.CPU_THREADS, exeflags=["--project=Loki.jl", "--heap-size-hint=4G"])
@everywhere using Loki

# Load in data
obs = from_fits(["data_ifualign/Level3_ch1-long_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch1-medium_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch1-short_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch2-long_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch2-medium_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch2-short_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch3-long_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch3-medium_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch3-short_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch4-long_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch4-medium_s3d_corr_nostripe.fits",
                 "data_ifualign/Level3_ch4-short_s3d_corr_nostripe.fits"],
                 0.597)

channel = 234
# channel = 34
nm = replace(obs.name, " " => "_") 
name = nm * "_ch$(channel)_psftemp_2comp_cubic_SCALECH_extmap_fullrun"

if isfile("$nm.channel$(channel).rest_frame.fits")
    obs = from_fits(["$nm.channel$(channel).rest_frame.fits"], obs.z)
else
    # Generate the PSF models
    generate_psf_model!(obs, "./psfs_stars")
    # Convert to rest-frame wavelength vector, and mask out bad spaxels
    correct!(obs)
    # Concatenate the subchannels of each channel so that we have one cube for each channel
    for i_channel in 2:4
        for i_band in (:A, :B, :C)
            splinefit_psf_model!(obs.channels[Symbol(i_band, i_channel)], 100)
        end
    end

    scale_factors = combine_channels!(obs, [:A2, :B2, :C2, :A3, :B3, :C3, :A4, :B4, :C4], out_id=channel, order=1, rescale_channels=nothing,
        adjust_wcs_headerinfo=true, match_psf=false, rescale_limits=(0.0, Inf), rescale_snr=10., extract_from_ap=0., max_λ=17.2,
        rescale_all_psf=false, scale_psf_only=true, fix_4c=true)

    # Fit cubic spline to the PSF
    # splinefit_psf_model!(obs.channels[channel], 100)
    # Rotate to sky axes
    rotate_to_sky_axes!(obs.channels[channel])
    # Interpolate NaNs in otherwise good spaxels
    interpolate_nans!(obs.channels[channel])
    calculate_statistical_errors!(obs.channels[channel], 20, 5, 3.0)
    # Save the pre-processed data as a FITS file so it can be quickly reloaded later
    save_fits(".", obs, channel)

    obs = from_fits(["$nm.channel$(channel).rest_frame.fits"], obs.z)
end

# Generate the 3D nuclear template using the PSF models
nuc_temp = generate_nuclear_template(obs.channels[channel], 0.)

plot_range = [(5.0, 5.1), (5.3, 5.65), (6.08, 6.13), (6.89, 7.0), (7.45, 8.05), (8.79, 9.19), (9.46, 9.86), (10.31, 10.71),
    (12.2, 12.35), (12.6, 13.0), (14.0, 14.6), (15.4, 15.7)]
linemask_overrides = [(5.04, 5.06), (5.32, 5.35), (5.435, 5.467), (5.489, 5.52), (5.59, 5.63), (6.098, 6.117), (6.895, 6.926), (6.96, 7.007),
    (7.43, 7.47), (7.62, 7.68), (8.01, 8.04), (8.955, 9.04), (9.492, 9.58), (9.61, 9.71), (10.47, 10.55),
    (12.257, 12.298), (12.77, 12.86), (14.256, 14.385), (15.49, 15.61), (16.95, 17.10)]
user_mask = []
line_test_lines = [
    ["H200_S1*", "H200_S2*", "H200_S3", "H200_S4*", "H200_S5*", "H200_S6*", "H200_S7*", "H200_S8*"], 
    ["ArII_6985*", "NeII_12813", "FeII_5340*"],
    ["SIV_10511", "ArIII_8991*", "NeIII_15555*"],
    ["NeVI_7652", "FeVIII_5447*", "MgVII_5503*", "MgV_5609*", "MgVII_9009*", "NeV_14322"],
]

fit_stellar_continuum = false
use_pah_templates = true
# use_pah_templates = false
custom_ext_template = readdlm("ext_template_phoenix.txt", '\t', Float64, '\n')

# extinction_map = nothing
# extinction map from channel 3
extinction_map = readdlm("PHOENIX.EXTINCTION.FROM_H2.BLUR1.CSV", ',', Float64, '\n')
extinction_map = Matrix(extinction_map')

# ch3param = FITS("GOOD6_output_Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir/Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir_parameter_maps.fits")
# ch3param = FITS("GOOD6_output_Phoenix_Cluster_BCG_ch3_NOpsftemp_2comp_SCALECH_mir/Phoenix_Cluster_BCG_ch3_NOpsftemp_2comp_SCALECH_mir_parameter_maps.fits")
# ch3param = FITS("GOOD7_output_Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir/Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir_parameter_maps.fits")
ch3param = FITS("Phoenix_Cluster_BCG.channel3.rest_frame.fits")

# reproject onto the current wcs
wcs3 = WCS.from_header(read_header(ch3param[2], String))[1]
wcs3 = WCSTransform(2; cdelt=wcs3.cdelt[1:2], crpix=wcs3.crpix[1:2], crval=wcs3.crval[1:2], ctype=wcs3.ctype[1:2], cunit=wcs3.cunit[1:2],
                    pc=wcs3.pc[1:2, 1:2], radesys=wcs3.radesys)
wcsnow = obs.channels[channel].wcs
wcsnow = WCSTransform(2; cdelt=wcsnow.cdelt[1:2], crpix=wcsnow.crpix[1:2], crval=wcsnow.crval[1:2], ctype=wcsnow.ctype[1:2], cunit=wcsnow.cunit[1:2],
                    pc=wcsnow.pc[1:2, 1:2], radesys=wcsnow.radesys)
shapenow = size(obs.channels[channel].I)[1:2]
extinction_map[.~isfinite.(extinction_map)] .= 0.
extinction_map, _ = reproject((extinction_map, wcs3), wcsnow, shape_out=shapenow, order=1)
extinction_map[.~isfinite.(extinction_map)] .= 0.

# Create the cube fitting object
# second iteration fit allowing extinction to vary freely
cube_fitter = CubeFitter(
    obs.channels[channel], 
    obs.z, 
    name; 
    parallel=true, 
    # parallel_strategy="pmap",
    plot_spaxels=:pyplot, 
    plot_maps=true, 
    save_fits=true,
    extinction_curve="custom", 
    custom_ext_template=custom_ext_template,
    extinction_F_test=true,
    extinction_screen=true, 
    use_pah_templates=use_pah_templates, 
    fit_sil_emission=false, 
    fit_stellar_continuum=fit_stellar_continuum, 
    fit_temp_multexp=false,
    linemask_delta=20, 
    linemask_n_inc_thresh=5, 
    linemask_thresh=3., 
    linemask_overrides=linemask_overrides, 
    save_full_model=true, 
    map_snr_thresh=3., 
    user_mask=user_mask, 
    fit_all_global=false,
    # templates=obs.channels[channel].psf_model,
    templates=nuc_temp, 
    template_names=["nuclear"], 
    extinction_map=extinction_map,
    subtract_cubic_spline=true,
    line_test_lines=line_test_lines, 
    line_test_threshold=0.003,
    lock_hot_dust=false,
    sort_line_components=:flux,
    # p_init_cube="output_Phoenix_Cluster_BCG_ch0_multichannel_234_blurred_psftemp_70pm10_mir",
    plot_range=plot_range
)

# Fit the cube
try
    fit_cube!(cube_fitter)
finally
    rmprocs(procs)
end
