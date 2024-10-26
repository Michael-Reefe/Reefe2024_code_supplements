using Distributed
using Pkg
using FITSIO
using Photometry
using DelimitedFiles
using Dierckx
using WCS
using Reproject
using NaNStatistics

procs = addprocs(Sys.CPU_THREADS, exeflags=["--project=Loki.jl" "--heap-size-hint=4G"])
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

channel = parse(Int, ARGS[1])
nm = replace(obs.name, " " => "_") 
name = nm * "_ch$(channel)_psftemp_2comp_W80FIX"

if isfile("$nm.channel$channel.rest_frame.fits")
    obs = from_fits(["$nm.channel$channel.rest_frame.fits"], obs.z)
else
    # Generate the PSF models
    generate_psf_model!(obs, "./psfs_stars")
    # generate_psf_model!(obs, "../Loki/src/templates/webbpsf", interpolate_leak_artifact=false)
    # Convert to rest-frame wavelength vector, and mask out bad spaxels
    correct!(obs)
    # Fit cubic spline to the PSF
    for band in (:A, :B, :C)
        chband = Symbol(band, channel)
        splinefit_psf_model!(obs.channels[chband], 100)
    end
    # Concatenate the subchannels of each channel so that we have one cube for each channel
    scale_factors = combine_channels!(obs, channel, out_id=channel, order=1, rescale_channels=nothing,
        adjust_wcs_headerinfo=true, rescale_limits=(0.0, Inf), rescale_snr=10., extract_from_ap=0., max_λ=17.2)

    # Rotate to sky axes
    rotate_to_sky_axes!(obs.channels[channel])
    # Interpolate NaNs in otherwise good spaxels
    interpolate_nans!(obs.channels[channel])
    calculate_statistical_errors!(obs.channels[channel], 20, 5, 3.0)
    # Save the pre-processed data as a FITS file so it can be quickly reloaded later
    save_fits(".", obs, channel)
    # Load it back in
    obs = from_fits(["$nm.channel$channel.rest_frame.fits"], obs.z)
end

# Generate the 3D nuclear template using the PSF models
nuc_temp = generate_nuclear_template(obs.channels[channel], 0.)

use_pah_templates = true
if channel == 3
    plot_range = [(7.45, 8.05), (8.79, 9.19), (9.46, 9.86), (10.31, 10.71)]
    linemask_overrides = [(7.43, 7.47), (7.61, 7.69), (8.01, 8.04), (8.955, 9.07), (9.492, 9.566), (9.61, 9.71), (10.472, 10.547)]
    line_test_lines = [["H200_S3", "H200_S4*"], ["SIV_10511", "ArIII_8991*"], ["NeVI_7652", "MgVII_9009*", "FeVII_9527*"]]
    # line_test_lines = []
    user_mask = []
    # line_test_lines = []
    # user_mask = [[9.54, 9.64], [9.69, 9.8]]
    # user_mask = [[9.6, 10.0]]
    fit_stellar_continuum = false
elseif channel == 2
    plot_range = [(5.0, 5.1), (5.3, 5.65), (6.08, 6.13), (6.89, 7.0)]
    linemask_overrides = [(5.04, 5.06), (5.32, 5.35), (5.435, 5.467), (5.489, 5.52), (5.59, 5.63), (6.098, 6.117), (6.895, 6.926), (6.96, 7.007)]
    user_mask = []
    line_test_lines = [["H200_S5", "H200_S6*", "H200_S7*", "H200_S8*"], ["FeVIII_5447*", "MgV_5609", "MgVII_5503*"], ["FeII_5340*", "ArII_6985"]]
    # line_test_lines = []
    fit_stellar_continuum = false
elseif channel == 4
    plot_range = [(12.2, 12.35), (12.6, 13.0), (14.0, 14.6), (15.4, 15.7), (16.8, 17.2)]
    linemask_overrides = [(12.257, 12.298), (12.77, 12.86), (14.256, 14.385), (15.49, 15.61), (16.95, 17.10)]
    user_mask = []
    line_test_lines = [["H200_S2", "H200_S1*"], ["NeII_12813", "NeIII_15555"], ["NeV_14322"]]
    # line_test_lines = []
    fit_stellar_continuum = false
elseif channel == 1
    plot_range = [(3.9, 4.1), (4.45, 4.6)]
    linemask_overrides = [(3.92, 3.95), (4.03, 4.07), (4.468, 4.551)]
    user_mask = []
    line_test_lines = [["MgIV_4867", "ArVI_4529"], ["SiIX_3936"]]
    # line_test_lines = []
    fit_stellar_continuum = true
    use_pah_templates = false
end

# extinction map from channel 3
extinction_map = readdlm("PHOENIX.EXTINCTION.FROM_H2.BLUR1.CSV", ',', Float64, '\n')
extinction_map = Matrix(extinction_map')

if channel != 3
    # ch3param = FITS("GOOD6_output_Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir/Phoenix_Cluster_BCG_ch3_psftemp_2comp_SCALECH_mir_parameter_maps.fits")
    # ch3param = FITS("GOOD6_output_Phoenix_Cluster_BCG_ch3_NOpsftemp_2comp_SCALECH_mir/Phoenix_Cluster_BCG_ch3_NOpsftemp_2comp_SCALECH_mir_parameter_maps.fits")
    # ch3param = FITS("output_Phoenix_Cluster_BCG_ch3_psftemp_2comp_W80FIX_mir/Phoenix_Cluster_BCG_ch3_psftemp_2comp_W80FIX_mir_parameter_maps.fits")
    ch3param = FITS("Phoenix_Cluster_BCG.channel3.rest_frame.fits")

    # extinction_map = read(ch3param["EXTINCTION_TAU_9_7"])
    # extinction_map = read(ch3param["EXTINCTION.TAU_9_7"])
    # extinction_map = cat(read(ch3param["EXTINCTION_N_OLI"]), read(ch3param["EXTINCTION_N_PYR"]), read(ch3param["EXTINCTION_N_FOR"]), dims=3)

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
    # for i in axes(extinction_map, 3)
        # e, _ = reproject((extinction_map[:,:,i], wcs3), wcsnow, shape_out=shapenow, order=1)
        # extinction_map[:,:,i] .= e
    # end
    extinction_map[.~isfinite.(extinction_map)] .= 0.

end
custom_ext_template = readdlm("ext_template_phoenix.txt", '\t', Float64, '\n')

# Create the cube fitting object
# second iteration fit allowing extinction to vary freely
cube_fitter = CubeFitter(
    obs.channels[channel], 
    obs.z, 
    name; 
    parallel=true, 
    parallel_strategy="pmap",
    plot_spaxels=:pyplot, 
    plot_maps=true, 
    save_fits=true,
    extinction_curve="custom",
    extinction_screen=true, 
    custom_ext_template=custom_ext_template,
    use_pah_templates=use_pah_templates, 
    fit_sil_emission=false, 
    fit_stellar_continuum=fit_stellar_continuum, 
    linemask_delta=20, 
    linemask_n_inc_thresh=5, 
    linemask_thresh=3., 
    linemask_overrides=linemask_overrides, 
    save_full_model=true, 
    map_snr_thresh=3., 
    user_mask=user_mask, 
    fit_all_global=true,
    templates=nuc_temp, 
    template_names=["nuclear"], 
    # tie_template_amps=true,
    extinction_map=extinction_map,
    subtract_cubic_spline=true,
    line_test_lines=line_test_lines, 
    line_test_threshold=0.003,  # require a 3-sigma second component (99.7% -> 1-0.997=0.003) 
    sort_line_components=:flux, 
    # p_init_cube="output_Phoenix_Cluster_BCG_ch0_multichannel_234_blurred_psftemp_70pm10_mir",
    plot_range=plot_range,
    # lock_hot_dust=true,
    # F_test_ext=true
)

# Fit the cube
try
    fit_cube!(cube_fitter)
finally
    rmprocs(procs)
end
