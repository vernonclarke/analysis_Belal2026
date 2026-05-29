# STANDARDIZED NWB CONVERSION FUNCTION
#' Convert ABF files to NWB format with standardized workflow
#' 
#' @param experimental_dict list containing experiment metadata
#' @param username system username for path construction
#' @param path_analysis base path for analysis files
#' @param identifier subfolder identifier (e.g., 'Figure 1', 'Figure 3/CRISPR')
#' @param default_stimulation default stimulation time if auto-detection fails (ms)
#' @param folder_name_map optional named list mapping experiment names to NWB folder names
#' @param get_metadata_fn function that returns metadata list for a given experiment_name
#' @param dandi_compliant if TRUE puts all experiments into one folder (identifier) else into subfolders (not compliant)
#' 
convert_abf_to_nwb <- function(experimental_dict, 
                                username,
                                path_analysis,
                                identifier,
                                default_stimulation = 150,
                                default_baseline = 100,
                                folder_name_map = NULL,
                                get_metadata_fn,
                                dandi_compliant = TRUE) { 
  
  # Construct paths
  file_path3 <- file.path('/Users', username, path_analysis, 'Raw ABF data', identifier)
  # Apply folder_name_map to NWB path if provided
  nwb_identifier <- identifier
  if (!is.null(folder_name_map)) {
    for (old_name in names(folder_name_map)) {
      if (grepl(old_name, identifier, fixed = TRUE)) {
        nwb_identifier <- gsub(old_name, folder_name_map[[old_name]], identifier, fixed = TRUE)
        break
      }
    }
  }
  nwb_base_path <- file.path('/Users', username, path_analysis, 'NWB data', nwb_identifier)
    
  # Create base NWB directory
  if (!dir.exists(nwb_base_path)) {
    dir.create(nwb_base_path, recursive = TRUE)
    cat("Created NWB base directory:", nwb_base_path, "\n")
  }
  
  # Get experiment IDs
  expt_id <- names(experimental_dict)
  
  # Process each experiment
  for (iii in seq_along(expt_id)) {
    
    experiment_name <- expt_id[iii]
    cat("\n=== Processing experiment:", experiment_name, "===\n")
    
    # Determine output path based on DANDI compliance
    if (dandi_compliant) {
      # DANDI-compliant: all subjects directly at root
      dandiset_path <- nwb_base_path
      cat("Using DANDI-compliant flat structure\n")
    } else {
      # Legacy: separate folders for each experiment
      if (!is.null(folder_name_map) && experiment_name %in% names(folder_name_map)) {
        nwb_folder_name <- folder_name_map[[experiment_name]]
      } else {
        nwb_folder_name <- experiment_name
      }
      dandiset_path <- file.path(nwb_base_path, nwb_folder_name)
      if (!dir.exists(dandiset_path)) {
        dir.create(dandiset_path, recursive = TRUE)
        cat("Created experiment directory:", dandiset_path, "\n")
      }
    }
    
    # Get folders for this dataset
    folders <- names(experimental_dict[[iii]])
    
    # Loop through all folders
    for (ii in seq_along(folders)) {
      
      folder <- folders[ii]
      cat("\n--- Processing folder:", folder, "---\n")
      
      # ABF folder
      # Check if identifier already contains the experiment name
      if (endsWith(identifier, experiment_name)) {
        abf_path <- file.path(file_path3, folder)
      } else {
        abf_path <- file.path(file_path3, experiment_name, folder)
      }
      
      if (!dir.exists(abf_path)) {
        cat("WARNING: ABF path does not exist:", abf_path, "\n")
        next
      }
      
      setwd(abf_path)
      
      # Find ABF files
      abf_files <- list.files(path = '.', pattern = '\\.[Aa][Bb][Ff]$', full.names = FALSE)
      
      if (length(abf_files) == 0) {
        cat("WARNING: No ABF files found in:", abf_path, "\n")
        next
      }
      
      cat("Found", length(abf_files), "ABF file(s)\n")
      
      # Read ABF data and convert
      tryCatch({
        abf_dataset <- readABFs(abf_files)
        abf_file <- abf_files[1]
        
        # Get timestamp
        session_start_time <- time_stamp(abf_dataset)
        cat("ABF recording time:", session_start_time, "\n")
        
        # Generate subject ID and filename
        base_id <- sub('\\.abf$', '', abf_file, ignore.case = TRUE)
        subject_id <- paste0('m', base_id)
        
        nwb_filename <- generate_dandi_nwb_filename(subject_id, session_id = NULL, modality = 'icephys')
        
        # Create DANDI-compliant subject folder
        subject_folder <- file.path(dandiset_path, paste0('sub-', subject_id))
        if (!dir.exists(subject_folder)) {
          dir.create(subject_folder, recursive = TRUE)
        }
        
        # Path for NWB output
        nwb_output_path <- file.path(subject_folder, nwb_filename)
        
        # Get metadata from experimental dict
        traces2average <- experimental_dict[[iii]][[folder]]$traces2average
        level_names <- experimental_dict[[iii]][[folder]]$levels
        traces2save <- as.integer(unlist(traces2average) - 1)
        
        # Calculate level_column_mapping
        level_column_mapping <- list()
        col_start <- 0
        for (i in seq_along(level_names)) {
          num_traces_in_level <- length(traces2average[[i]])
          if (num_traces_in_level > 0) {
            level_column_mapping[[level_names[i]]] <- seq(col_start, col_start + num_traces_in_level - 1)
            col_start <- col_start + num_traces_in_level
          }
        }
        
        cat("Level column mapping for", folder, ":\n")
        print(level_column_mapping)
        
        # Calculate dt (sampling interval in ms)
        dt <- abf_dataset$samplingIntervalInSec * 1000

        # Calculate recording duration
        recording_duration_sec <- (max(abf_dataset$header$sweepStartInPts) + abf_dataset$header$sweepLengthInPts) * dt / 1000

        # Set baseline
        baseline <- default_baseline  # ms
        
        # Auto-detect stimulation time from TTL channel
        ttl_channel <- 3
        stim_result <- detect_stimulation_from_ttl(abf_dataset, ttl_channel, dt)
        
        if (stim_result$success) {
          stimulation_time <- stim_result$time
          cat("Auto-detected stimulation time:", round(stimulation_time, 2), "ms\n")
        } else {
          stimulation_time <- default_stimulation
          cat("Using default stimulation time:", stimulation_time, "ms\n")
        }
        
        # Get experiment-specific metadata using provided function
        metadata <- get_metadata_fn(experiment_name)
        
        # Generate NWB file
        cat("Generating NWB file:", nwb_filename, "\n")
        
        abf2nwb(
          abf_dataset = abf_dataset,
          nwb_dataset = nwb_output_path,
          subject_id = subject_id,
          session_start_time = session_start_time,
          recording_duration = recording_duration_sec,
          species = metadata$species,
          age = metadata$age,
          sex = metadata$sex,
          genotype = metadata$genotype,
          cross = metadata$cross,
          virus = metadata$virus,
          virus_injection_site = metadata$virus_injection_site,
          experimenter = metadata$experimenter,
          institution = metadata$institution,
          lab = metadata$lab,
          session_description = metadata$session_description,
          experiment_description = metadata$experiment_description,
          keywords = metadata$keywords,
          device_description = metadata$device_description,
          traces2save = traces2save,
          levs = level_column_mapping,
          baseline = baseline,
          stimulation_time = stimulation_time,
          notes = metadata$notes
        )
        
        cat("Saved:", nwb_output_path, "\n")
        
      }, error = function(e) {
        cat("ERROR processing", folder, ":", conditionMessage(e), "\n")
      })
    }
  }
  
  cat("\n=== BATCH CONVERSION COMPLETE ===\n")
  cat("NWB files saved to:", nwb_base_path, "\n")
  
  invisible(nwb_base_path)
}

#' Print NWB File Metadata
#' 
#' Reads and displays all metadata from an NWB file including session info,
#' subject details, acquisition data, and processing modules.
#'
#' @param nwb_filepath Full path to the NWB file
#' @return prints metadata to console.
#' @examples
#' print_nwb_metadata("/path/to/file.nwb")
#' 
print_nwb_metadata <- function(nwb_filepath) {
  
  if (!file.exists(nwb_filepath)) {
    stop("NWB file not found: ", nwb_filepath)
  }
  
  # Pass filepath to Python
  py$nwb_filepath <- nwb_filepath
  
  # Read and print all metadata
  py_run_string("
from pynwb import NWBHDF5IO

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()

# Print all metadata
print('\\n=== SESSION INFO ===')
print(f'Session description: {nwb.session_description}')
print(f'Session start time: {nwb.session_start_time}')
print(f'Experimenter: {\", \".join(nwb.experimenter)}')
print(f'Institution: {nwb.institution}')
print(f'Lab: {nwb.lab}')

print('\\n=== SUBJECT INFO ===')
print(f'Subject ID: {nwb.subject.subject_id}')
print(f'Species: {nwb.subject.species}')
print(f'Age: {nwb.subject.age}')
print(f'Sex: {nwb.subject.sex}')
print(f'Genotype: {nwb.subject.genotype}')
print(f'Description: {nwb.subject.description}')

print('\\n=== EXPERIMENT INFO ===')
print(f'Experiment description: {nwb.experiment_description}')
print(f'Keywords: {nwb.keywords[:]}')

print('\\n=== DEVICES ===')
for device_name, device in nwb.devices.items():
    print(f'{device_name}:')
    print(f'  Description: {device.description}')

print('\\n=== ELECTRODES ===')
for elec_name, elec in nwb.icephys_electrodes.items():
    print(f'{elec_name}:')
    print(f'  Cell ID: {elec.cell_id}')
    print(f'  Description: {elec.description}')
    print(f'  Location: {elec.location}')
    print(f'  Resistance: {elec.resistance}')
    print(f'  Filtering: {elec.filtering}')

print('\\n=== ACQUISITION DATA ===')
print(f'Number of sweeps: {len(nwb.acquisition)}')
for name, series in nwb.acquisition.items():
    print(f'  Sweep {series.sweep_number}: {series.description}')
    print(f'    Rate: {series.rate:.1f} Hz')
    print(f'    Gain: {series.gain:.6f}')
    print(f'    Unit: {series.unit}')

print('\\n=== PROCESSING MODULES ===')
for module_name in nwb.processing.keys():
    print(f'Module: {module_name}')
    module = nwb.processing[module_name]
    for container_name, container in module.data_interfaces.items():
        if hasattr(container, 'data'):
            data_value = container.data[:]
            if container_name == 'level_column_mapping':
                import json
                json_str = ''.join([chr(int(c)) for c in data_value])
                data_value = json.loads(json_str)
            print(f'  {container_name}: {data_value}')

print('\\n=== FILE INFO ===')
print(f'Identifier: {nwb.identifier}')
print(f'File created: {nwb.file_create_date[0].strftime(\"%Y-%m-%d %H:%M:%S %Z\")}')

io.close()
", convert = FALSE)
  
  invisible(NULL)
}

# print_nwb_metadata <- function(nwb_filepath) {
  
#   if (!file.exists(nwb_filepath)) {
#     stop("NWB file not found: ", nwb_filepath)
#   }
  
#   # Pass filepath to Python
#   py$nwb_filepath <- nwb_filepath
  
#   # Read and print all metadata
#   py_run_string("
# from pynwb import NWBHDF5IO

# io = NWBHDF5IO(nwb_filepath, 'r')
# nwb = io.read()

# # Print all metadata
# print('\\n=== SESSION INFO ===')
# print(f'Session description: {nwb.session_description}')
# print(f'Session start time: {nwb.session_start_time}')
# print(f'Experimenter: {\", \".join(nwb.experimenter)}')
# print(f'Institution: {nwb.institution}')
# print(f'Lab: {nwb.lab}')

# print('\\n=== SUBJECT INFO ===')
# print(f'Subject ID: {nwb.subject.subject_id}')
# print(f'Species: {nwb.subject.species}')
# print(f'Age: {nwb.subject.age}')
# print(f'Sex: {nwb.subject.sex}')
# print(f'Genotype: {nwb.subject.genotype}')

# print('\\n=== CUSTOM METADATA ===')
# for key, value in nwb.fields.items():
#     if key not in ['session_description', 'session_start_time', 'experimenter', 'institution', 'lab', 'subject', 'processing', 'keywords', 'acquisition']:
#         if key == 'file_create_date':
#             formatted_date = value[0].strftime('%Y-%m-%d %H:%M:%S%z')
#             print(f'{key}: {formatted_date}')
#         else:
#             print(f'{key}: {str(value)}')

# print(f'\\nKeywords: {nwb.keywords[:]}')

# print('\\n=== ACQUISITION DATA ===')
# print(f'Number of sweeps: {len(nwb.acquisition)}')
# for name, series in nwb.acquisition.items():
#     print(f'  Sweep {series.sweep_number}: {series.description}')
#     print(f'    Rate: {series.rate:.1f} Hz, Gain: {series.gain:.6f}, Unit: {series.unit}')

# print('\\n=== PROCESSING MODULES ===')
# for module_name in nwb.processing.keys():
#     print(f'Module: {module_name}')
#     module = nwb.processing[module_name]
#     for container_name, container in module.data_interfaces.items():
#         if hasattr(container, 'data'):
#             data_value = container.data[:]
#             if container_name == 'level_column_mapping':
#                 import json
#                 json_str = ''.join([chr(int(c)) for c in data_value])
#                 data_value = json.loads(json_str)
#             print(f'  - {container_name}: {data_value}')
#         else:
#             print(f'  - {container_name}: {str(container)}')

# io.close()
# ", convert = FALSE)
  
#   invisible(NULL)
# }


load_nwb_averages <- function(nwb_filepath) {
  
  py$nwb_filepath <- nwb_filepath
  
  py_run_string("
from pynwb import NWBHDF5IO
import numpy as np
import json

io = NWBHDF5IO(nwb_filepath, 'r')
nwb = io.read()

baseline = float(nwb.processing['icephys']['baseline_duration'].data[0])
stim_time = float(nwb.processing['icephys']['stimulation_time'].data[0])
samplingIntervalInSec = float(nwb.processing['icephys']['samplingIntervalInSec'].data[0])

# Decode level_column_mapping
level_mapping_raw = nwb.processing['icephys']['level_column_mapping'].data[:]
level_mapping_json = ''.join([chr(int(c)) for c in level_mapping_raw])
level_column_mapping = json.loads(level_mapping_json)

level_names = list(level_column_mapping.keys())

# Get acquisition keys
acquisition_keys = sorted(nwb.acquisition.keys(), key=lambda x: int(x))

# Get rate from first acquisition series
first_key = acquisition_keys[0]
rate = float(nwb.acquisition[first_key].rate)

# Extract all traces and convert from amperes to pA
traces = []
for key in acquisition_keys:
    trace_data = np.array(nwb.acquisition[key].data[:])
    trace_pA = trace_data * 1e12  # amperes to pA
    traces.append(trace_pA)

io.close()
")
  
  # extract variables from Python
  baseline <- py$baseline
  stimulation_time <- py$stim_time
  level_column_mapping <- py$level_column_mapping
  level_names <- unlist(py$level_names)
  rate <- py$rate
  
  # get traces as a matrix
  traces_list <- py$traces
  traces <- do.call(cbind, traces_list)
  
  # calculate indices
  nlevels <- length(level_names)
  dt <- py$samplingIntervalInSec * 1000
  idx1 <- (stimulation_time - baseline) / dt  
  idx2 <- baseline / dt
  
  # calculate averaged traces for each level
  trace_averages <- sapply(1:nlevels, function(ii) {
    level_traces <- traces[, level_column_mapping[[level_names[ii]]] + 1, drop = FALSE]
    traces_subset <- level_traces[idx1:nrow(level_traces), , drop = FALSE]
    baselines <- colMeans(traces_subset[1:idx2, , drop = FALSE])
    baseline_corrected <- sweep(traces_subset, 2, baselines)
    rowMeans(baseline_corrected)
  })
  
  colnames(trace_averages) <- level_names
  
  # return as data.frame (analogous to out1$`single examples`[-1])
  as.data.frame(trace_averages)
}


indices2traces <- function(indices, nmax=5) {
  indices <- sort(unique(indices))
  n_indices <- length(indices)
  if (n_indices == 0) return(list())
  # Calculate non-overlapping group sizes (distance from previous index)
  possible_sizes <- sapply(seq_along(indices), function(i) {
    if (i == 1) {
      max(1, indices[1] - 1)
    } else {
      max(1, indices[i] - indices[i-1] - 1)
    }
  })
  group_size <- min(nmax, possible_sizes)
  if (group_size < 1)
    stop("Not enough room before one or more indices for at least 1 trace.")

  out <- vector("list", n_indices)
  prev_max <- 0
  for (i in seq_along(indices)) {
    idx <- indices[i]
    if (i == 1) {
      group <- seq(idx - group_size, idx - 1)
    } else {
      # Prevent overlap with previous group/index
      start <- max(idx - group_size, prev_max + 1, indices[i-1] + 1)
      group <- seq(start, idx - 1)
      if (length(group) > group_size) group <- tail(group, group_size)
    }
    group <- group[group > 0]
    out[[i]] <- group
    if (length(group)) prev_max <- max(group)
  }
  out
}

metadata_extractor <- function(abf_dataset, electrode_resistance='3-5 MOhm', channel_index=1, traces2save=NULL, levs=NULL) {
  
  # extract all header information from the abf file
  metadata <- extract_metadata (abf_dataset)

  # extract useful information
  unit <- metadata$channelUnits[channel_index]
  gain_headstage <- metadata$header$fInstrumentScaleFactor[channel_index]
  gain_additional <- metadata$header$fTelegraphAdditGain[channel_index]
  gain_total_mV_per_unit <- gain_headstage * gain_additional * 1000
  filtering <- metadata$sections$ADCsec[[channel_index]]$fSignalLowpassFilter
  resistance_comp <- metadata$sections$ADCsec[[channel_index]]$fTelegraphAccessResistance
  capacitance_comp <- metadata$sections$ADCsec[[channel_index]]$fTelegraphMembraneCap
  rate <- 1 / metadata$samplingIntervalInSec

  sweep_starts <- round(metadata$header$sweepStartInPts)
  sweep_length <- round(metadata$header$sweepLengthInPts)
  num_sweeps <- round(metadata$header$lActualEpisodes)

  # if traces2save is NULL all traces are converted 
  # 
  if (!is.null(traces2save)) {
    traces2save <- as.integer(traces2save - 1)
    sweep_starts <- sweep_starts[traces2save + 1]
    num_sweeps <- length(traces2save)
  }

  trace_data <- abf_dataset$data
  if (!is.null(trace_data[[1]]) && is.matrix(trace_data[[1]])) {
    all_traces <- lapply(trace_data, function(mat) mat[, channel_index])
  } else {
    stop('Invalid or missing trace data')
  }

  return(list(
    data = all_traces,
    unit = unit,
    gain_headstage = gain_headstage,
    gain_additional = gain_additional,
    gain_total_mV_per_unit = gain_total_mV_per_unit,
    filtering = filtering,
    electrode_resistance = electrode_resistance,
    resistance_comp = resistance_comp,
    capacitance_comp = capacitance_comp,
    rate = rate,
    sweep_starts = sweep_starts,
    sweep_length = sweep_length,
    num_sweeps = num_sweeps,
    traces2save = traces2save,
    levs = levs
  ))
}

generate_dandi_nwb_filename <- function(subject_id,
                                        session_id = NULL,
                                        acquisition = NULL,
                                        task = NULL,
                                        run = NULL,
                                        suffix = 'ephys',
                                        extension = '.nwb') {
  if (missing(subject_id) || !nzchar(subject_id)) {
    stop('A non-empty subject_id is required (e.g., "Mouse001").')
  }

  parts <- paste0('sub-', subject_id)

  if (!is.null(session_id)) {
    parts <- paste(parts, paste0('ses-', session_id), sep = '_')
  }

  if (!is.null(acquisition)) {
    parts <- paste(parts, paste0('acq-', acquisition), sep = '_')
  }

  if (!is.null(task)) {
    parts <- paste(parts, paste0('task-', task), sep = '_')
  }

  if (!is.null(run)) {
    parts <- paste(parts, paste0('run-', run), sep = '_')
  }

  filename <- paste0(parts, '_', suffix, extension)
  return(filename)
}

# abf2nwb <- function(abf_dataset, nwb_dataset,
#                     subject_id = 'Mouse001',
#                     species = 'Mus musculus',
#                     age = 'P30',
#                     sex = 'M',
#                     genotype = 'WT',
#                     cross = NULL,
#                     virus = NULL,
#                     virus_injection_site = NULL,
#                     experimenter = 'YourName',
#                     institution = 'YourInstitution',
#                     lab = 'YourLab',
#                     session_description = 'Whole-cell voltage clamp',
#                     device_name = 'Multiclamp 700B Patch clamp amplifier',
#                     electrode_resistance = '3-5 MOhm',
#                     location = 'Striatum',
#                     channel_index = 1,
#                     traces2save = NULL,
#                     levs = NULL) {

#   py_data <- metadata_extractor(
#     abf_dataset=abf_dataset,
#     electrode_resistance=electrode_resistance,
#     channel_index=channel_index,
#     traces2save=traces2save,
#     levs=levs
#   )

#   reticulate::py_run_string('
# import numpy as np
# from pynwb import NWBFile, NWBHDF5IO, TimeSeries, ProcessingModule
# from pynwb.file import Subject
# from pynwb.icephys import IntracellularElectrode, VoltageClampSeries
# from datetime import datetime

# def create_dandi_nwb(data_list, nwb_dataset,
#                      subject_id, species, age, sex, genotype,
#                      cross, virus, virus_injection_site,
#                      experimenter, institution, lab, session_description,
#                      device_name, location, traces2save):

#     rate = float(data_list["rate"])
#     unit = data_list["unit"]
#     gain_total = float(data_list["gain_total_mV_per_unit"])
#     gain_headstage = float(data_list["gain_headstage"])
#     gain_additional = float(data_list["gain_additional"])
#     filtering = str(data_list["filtering"]) + " Hz low-pass"
#     resistance = data_list["electrode_resistance"]
#     all_data = data_list["data"]

#     sweep_starts = np.array(data_list["sweep_starts"], dtype=np.int64)
#     sweep_length = np.array([data_list["sweep_length"]], dtype=np.int64)
#     num_sweeps = np.array([data_list["num_sweeps"]], dtype=np.int64)

#     n_traces = len(all_data)
#     traces = range(n_traces) if traces2save is None else traces2save

#     desc_parts = []
#     if cross:
#         desc_parts.append(f"Cross: {cross}")
#     if virus:
#         desc_parts.append(f"Virus: {virus}")
#     if virus_injection_site:
#         desc_parts.append(f"Injection site: {virus_injection_site}")
#     desc = " | ".join(desc_parts) if desc_parts else None

#     nwbfile = NWBFile(
#         session_description=session_description,
#         identifier="NWB_" + subject_id,
#         session_start_time=datetime.now(),
#         experimenter=[experimenter],
#         institution=institution,
#         lab=lab
#     )

#     subject = Subject(
#         subject_id=subject_id,
#         species=species,
#         age=age,
#         sex=sex,
#         genotype=genotype,
#         description=desc
#     )
#     nwbfile.subject = subject

#     device = nwbfile.create_device(name=device_name)

#     electrode = IntracellularElectrode(
#         name="elec0",
#         device=device,
#         description="Whole-cell patch-clamp electrode | Unit: {} | Headstage gain: {:.5g} V/{} | Additional gain: {:.2f}× | Total gain: {:.2f} mV/{}".format(
#             unit, gain_headstage, unit, gain_additional, gain_total, unit),
#         filtering=filtering,
#         location=location,
#         resistance=resistance
#     )
#     nwbfile.add_icephys_electrode(electrode)

#     traces = range(n_traces) if traces2save is None else traces2save

#     for j, i in enumerate(traces):
#         trace = all_data[i]
#         series = VoltageClampSeries(
#             name="response_{}".format(j),
#             data=trace,
#             rate=rate,
#             starting_time=0.0,
#             electrode=electrode,
#             gain=gain_total,
#             capacitance_slow=data_list["capacitance_comp"],
#             resistance_comp_correction=data_list["resistance_comp"],
#             stimulus_description="",
#             sweep_number=int(i),
#             unit=unit
#         )
#         nwbfile.add_acquisition(series)

#     # Add sweep info as a processing module
#     sweep_module = ProcessingModule(name="icephys", description="Sweep metadata from ABF header")
#     sweep_module.add_data_interface(TimeSeries(name="sweepStartInPts", data=sweep_starts, unit="pts", rate=1.0, starting_time=0.0))
#     sweep_module.add_data_interface(TimeSeries(name="sweepLengthInPts", data=sweep_length, unit="pts", rate=1.0, starting_time=0.0))
#     sweep_module.add_data_interface(TimeSeries(name="lActualEpisodes", data=num_sweeps, unit="count", rate=1.0, starting_time=0.0))

#     # Add levs and traces2save as additional metadata
#     if "levs" in data_list and data_list["levs"] is not None:
#         levs_arr = np.array(data_list["levs"], dtype=np.int64)
#         levs_series = TimeSeries(
#             name="levs",
#             data=levs_arr,
#             unit="none",
#             rate=1.0,
#             starting_time=0.0,
#             description="Levels for each trace"
#         )
#         sweep_module.add_data_interface(levs_series)

#     if "traces2save" in data_list and data_list["traces2save"] is not None:
#         traces_arr = np.array(data_list["traces2save"], dtype=np.int64)
#         traces_series = TimeSeries(
#             name="traces2save",
#             data=traces_arr,
#             unit="none",
#             rate=1.0,
#             starting_time=0.0,
#             description="Indices of traces to save"
#         )
#         sweep_module.add_data_interface(traces_series)

#     nwbfile.add_processing_module(sweep_module)

#     with NWBHDF5IO(nwb_dataset, "w") as io:
#         io.write(nwbfile)
# ')

#   reticulate::py$create_dandi_nwb(
#     py_data, nwb_dataset,
#     subject_id, species, age, sex, genotype,
#     cross, virus, virus_injection_site,
#     experimenter, institution, lab, session_description,
#     device_name, location, py_data$traces2save
#   )
# }

ABF_batch_analysis <- function(experimental_dict, file_path2, file_path3, baseline=100, stimulation=NULL, 
  graph_settings=list(width=6, height=8, xlim=NA, xlim1=NA, ylim1=NA, ylim2=NA, xlab='time (minutes)', 
    ylab='|PSC| (pA)', colors='#CD5C5C', xmajor_tick=10, ymajor_tick=100, cex.points=0.8),
  protocol=c('PSC'=1, 'hp'=2, 'TTL'=3), silent=FALSE, force_batch=TRUE, save=FALSE) {
  
  expt_id <- names(experimental_dict)
  ABF_summary <- lapply(1:length(expt_id), function(iii) {
    
    file_path5 <- paste0(file_path3, '/', expt_id[iii])
    setwd(file_path5)

    file_path6 <- sub('/Raw ABF data/', '/Raw ABF data summaries/', file_path3)
    if (save && !dir.exists(file_path6)) {
          dir.create(file_path6, recursive = TRUE)
      }

    new_path <- paste0(file_path6 , '/', expt_id[iii])
      if (save && !dir.exists(new_path)) {
        dir.create(new_path, recursive = TRUE)
      }

    # folders <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)
    folders <- names(experimental_dict[[expt_id[[iii]]]])
    
    if (length(folders) == 0) {
      out <- list()
    } else {
      out <- lapply(seq_along(folders), function(ii) {
        folder <- folders[ii]
        setwd(folder)

        xlsx_path1 <- paste0(new_path, '/', 'xlsx')
        if (save && !dir.exists(xlsx_path1)) {
          dir.create(xlsx_path1, recursive = TRUE)
        }
        xlsx_path2 <- paste0(xlsx_path1, '/', folder)
        if (save && !dir.exists(xlsx_path2)) {
          dir.create(xlsx_path2, recursive = TRUE)
        }

        svg_path1 <- paste0(new_path, '/', 'svg')
        if (save && !dir.exists(svg_path1)) {
          dir.create(svg_path1, recursive = TRUE)
        }
        svg_path2 <- paste0(svg_path1, '/', folder)
        if (save && !dir.exists(svg_path2)) {
          dir.create(svg_path2, recursive = TRUE)
        }

        entry <- experimental_dict[[expt_id[[iii]]]][[folder]]
        
        if ('traces2average' %in% names(entry) && !is.null(entry$traces2average)) {
          traces2average <- entry$traces2average
        }else if ('indices' %in% names(entry) && !is.null(entry$indices)) {
          traces2average <- indices2traces(indices=entry$indices,  nmax=5)
        }else{
          traces2average <- NULL
        }
        
        # abf_file=paste0(folders[ii], '.abf')
        abf_files <- list.files(path='.', pattern = '\\.[Aa][Bb][Ff]$', full.names=FALSE)

        message("Analysing ", paste0(folders[ii], '.abf'))
        tmp <- analyse_abf(abf_files, traces2average=traces2average, drug_application=entry$indices[-length(entry$indices)], 
          baseline=baseline, stimulation=stimulation, graph_settings=graph_settings, protocol=protocol, levels=entry$levels, file_path2=file_path2,
          xlsx_path=xlsx_path2, svg_path=svg_path2, silent=silent, force_batch=force_batch, save=save)

        # if (!force_batch){
        #   traces2average_final <- split(tmp$traces2average, tmp$traces2average[, 'level'])
        #   traces2average_final <- lapply(traces2average_final, function(x) x[x[, 'include'] == 1, 'trace'])
        #   tmp$'final traces2average' <- traces2average_final
        # }

        if (!force_batch) {
          trace_sel <- tmp[['trace selection']]

          if (!is.null(trace_sel) &&
              is.data.frame(trace_sel) &&
              all(c('Level', 'Accepted') %in% names(trace_sel))) {

            level_order <- entry$levels

            traces2average_final <- lapply(level_order, function(lvl) {
              rows <- which(trace_sel$Level == lvl)
              if (length(rows) == 0) return(integer(0))

              vals <- unlist(lapply(trace_sel$Accepted[rows], parse_trace_string), use.names = FALSE)
              vals <- as.integer(vals[!is.na(vals)])
              unique(vals)
            })

            names(traces2average_final) <- level_order
          } else {
            traces2average_final <- list()
          }

          tmp$`final traces2average` <- traces2average_final
        }

        setwd(file_path5) 
        tmp

      })

      names(out) <- folders
    }
    out
  })

  names(ABF_summary) <- expt_id
  setwd(file_path3)

  return(ABF_summary)
}


analyse_abf <- function(abf_files, traces2average, drug_application=NULL, baseline=100, stimulation=NULL,  
  protocol=c('PSC'=1, 'hp'=2, 'TTL'=3), levels = c('control', 'GABAzine'), file_path2 = '', xlsx_path = '.', svg_path = '.',
  graph_settings=list(width=6, height=8, xlim=NA, xlim1=NA, ylim1=NA, ylim2=NA, xlab='time (minutes)', ylab='|PSC| (pA)', colors='#CD5C5C', 
    xmajor_tick=10, ymajor_tick=100, cex.points=0.8), silent=FALSE, force_batch=FALSE, save=FALSE) {
       
  # abf_dataset <- readABF(abf_file)
  abf_files <- list.files(path='.', pattern = '\\.[Aa][Bb][Ff]$', full.names=FALSE)
  abf_dataset <-  readABFs(abf_files)
  abf_file <- abf_files[1] 
  abf_name <- tools::file_path_sans_ext(basename(abf_dataset$path))

  xlim <- graph_settings$xlim
  xlim1 <- graph_settings$xlim1
  ylim1 <- graph_settings$ylim1
  ylim2 <- graph_settings$ylim2
  xlab <- graph_settings$xlab
  ylab <- graph_settings$ylab
  width <- graph_settings$width
  height <- graph_settings$height
  colors <- graph_settings$colors
  xmajor_tick <- graph_settings$xmajor_tick
  ymajor_tick <- graph_settings$ymajor_tick
  cex.points <- graph_settings$cex.points

  dt <- abf_dataset$samplingIntervalInSec * 1000
  time= 0:dt:(length(abf_dataset$data[[1]][,1])-1)*dt

  I_data <- sapply(1:length(abf_dataset$data), function(ii) abf_dataset$data[[ii]][,protocol[['PSC']]])
  colnames(I_data) <- seq(dim(I_data)[2])
  rownames(I_data) <- seq(dim(I_data)[1])

  holding_potential <- sapply(1:length(abf_dataset$data), function(ii) abf_dataset$data[[ii]][,protocol[['hp']]])
  colnames(holding_potential) <- seq(dim(holding_potential)[2])
  rownames(holding_potential) <- seq(dim(holding_potential)[1])

  if ('TTL' %in% names(protocol) && is.null(stimulation)){
    TTL_pulse <- sapply(1:length(abf_dataset$data), function(ii) abf_dataset$data[[ii]][,protocol[['TTL']]])
    colnames(TTL_pulse) <- seq(dim(TTL_pulse)[2])
    rownames(TTL_pulse) <- seq(dim(TTL_pulse)[1])

    threshold <- 0.5
    pulse_on <- unname(which(TTL_pulse[,1] > threshold))

    if (length(pulse_on) > 0) {
      idx2 <- pulse_on[1]
      idx3 <- pulse_on[length(pulse_on)]
    } else {
      idx2 <- NA
      idx3 <- NA
    }

    stimulation_time <- idx2 * dt - dt

    out2 <- lapply(1:dim(TTL_pulse)[2], function(ii) peak.fun2(TTL_pulse[, ii], dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=5) )
    TTL_pulse2 <- sapply(1:length(out2), function(ii) out2[[ii]]$response)
  
  }else{
    if (is.null(stimulation)) {
      stop("With no TTL pulse found, 'stimulation' time must be specified")
    }else{
      stimulation_time <- stimulation
    }
  }

  out <- lapply(1:dim(I_data)[2], function(ii) peak.fun2(I_data[, ii], dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=5) )
  Apeak <- sapply(1:length(out), function(ii) out[[ii]]$peak)
  charge <- sapply(1:length(out), function(ii) out[[ii]]$charge)
  I_data2 <- sapply(1:length(out), function(ii) out[[ii]]$response)
  rownames(I_data2) <- seq(dim(I_data2)[1])
  colnames(I_data2) <- seq(dim(I_data2)[2])
  time <- seq(dim(I_data2)[1]) * dt  - dt
  holding_current <- sapply(1:length(out), function(ii) out[[ii]]$baseline)

  out1 <- lapply(1:dim(holding_potential)[2], function(ii) peak.fun2(holding_potential[, ii], dt=dt, stimulation_time=stimulation_time, baseline=baseline, smooth=5) )
  holding_potential <- sapply(1:length(out1), function(ii) out1[[ii]]$baseline)

  summary <- cbind(
    'time (minutes)' = 1:length(Apeak),
    'holding potential (mV)' = holding_potential,
    'holding current (pA)' = holding_current,
    'peak amplitude (pA)'  = -Apeak,
    'charge transfer (pC)'  = -charge
  )
  rownames(summary) <- 1:dim(summary)[1]

  Y <- list(summary[, 'peak amplitude (pA)'])
  X <- list(summary[, 'time (minutes)'])

  tmax <- max(summary[, 'time (minutes)'])
  xlim <- if (any(is.na(xlim))) c(0, 5 * ceiling(tmax / 5)) else xlim

  if (silent){
    
    split_include <- lapply(seq_along(traces2average), function(iii){
      select_traces(all_traces=I_data[,traces2average[[iii]], drop=FALSE], dt=dt, 
        amp=Apeak[traces2average[[iii]]], hc=holding_current[traces2average[[iii]]], 
        height=height, width=width, force_batch=force_batch)
    })

    result <- lapply(seq_along(split_include), function(iii) {
      mask <- split_include[[iii]]
      cols <- seq(mask)
      I_data2[,traces2average[[iii]], drop=FALSE][, cols[mask == 1], drop = FALSE]
    })
    average_traces <- sapply(result, function(mat) rowMeans(mat))
    single_examples <- cbind(time, average_traces)
    colnames(single_examples) <- c('time', levels)
    single_examples <- as.data.frame(single_examples)
    rownames(single_examples) <- seq_len(nrow(single_examples))
  
  }else{
    repeat {
      # select traces to include
      if (!force_batch) dev.new(width=width, height=height, noRStudioGD=TRUE)
      split_include <- lapply(seq_along(traces2average), function(iii){
        select_traces(all_traces=I_data[,traces2average[[iii]], drop=FALSE], dt=dt, 
          amp=Apeak[traces2average[[iii]]], hc=holding_current[traces2average[[iii]]], 
          height=height, width=width, force_batch=force_batch)
      })
      if (!force_batch) dev.off()

      result <- lapply(seq_along(split_include), function(iii) {
        mask <- split_include[[iii]]
        cols <- seq(mask)
        I_data2[,traces2average[[iii]], drop=FALSE][, cols[mask == 1], drop = FALSE]
      })
      average_traces <- sapply(result, function(mat) rowMeans(mat))
      single_examples <- cbind(time, average_traces)
      colnames(single_examples) <- c('time', levels)
      single_examples <- as.data.frame(single_examples)
      rownames(single_examples) <- seq_len(nrow(single_examples))
      
      if (any(is.na(ylim1))) {
        max_y <- ceiling(max(summary[, 'peak amplitude (pA)'], na.rm=TRUE) / 100) * 100
        ylim1 <- c(0, max_y)
      }

      for (lvl in rev(levels)) {     
        x <- single_examples[['time']]
        if (any(is.na(xlim1))) {
          xlim1 <- c(100*floor(min(x)/100), 100*ceiling(max(x)/100))
        }
        single_egs(x=x, y=single_examples[[lvl]], ylim=c(-ylim1[2],20), xlim=xlim1, 
          color='darkgrey', width=width, height=height)
      }

      # main plot
      dev.new(width=width, height=height, noRStudioGD=TRUE)
      layout(matrix(1:3, ncol = 1), heights = c(3, 2, 2))

      par(mar = c(3.5, 6.5, 2, 1), mgp = c(2.5, 0.5, 0))

      Y <- list(summary[, 'peak amplitude (pA)'])
      plot3(X, Y, xlim=xlim, ylim=ylim1, xlab='', ylab=ylab, cex.points=cex.points, 
            title= as.character(abf_name), colors=colors, ds=1, xmajor_tick=xmajor_tick, 
            show_xaxes=FALSE, ymajor_tick=ymajor_tick, show_spline=FALSE)

      if (length(traces2average) > 0) {
        for (iii in seq_along(split_include)) {
          mask <- split_include[[iii]]
          cols <- seq(mask)
          exY <- summary[traces2average[[iii]], 'peak amplitude (pA)'][cols[mask == 1], drop = FALSE]
          exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
          points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
        }
      }

      if (!is.null(drug_application)) {
        for (i in seq_along(drug_application)) {
          drug_box(
            box_times = c(drug_application[i], tmax),
            box_col = if (i %% 2 == 1) 'darkgrey' else 'lightgrey',
            box_border_col = NA,
            rel_box_height = 0.025,
            rel_y = 1 - (i - 1) * 0.1
          )
        }
      }

      Y <- list(summary[, 'holding current (pA)'])
      if (any(is.na(graph_settings$ylim2))) {
        min_y <- floor(min(summary[, 'holding current (pA)'], na.rm=TRUE) / 100) * 100
        ylim2 <- c(min_y, 0)
      }
      plot3(X, Y, xlim=xlim, ylim=ylim2, xlab='', ylab='holding current (pA)', cex.points=cex.points, 
            title='', colors=colors, ds=1, xmajor_tick=xmajor_tick, 
            ymajor_tick=ymajor_tick, show_xaxes=FALSE, show_spline=FALSE)

      if (length(traces2average) > 0) {
        for (iii in seq_along(split_include)) {
          mask <- split_include[[iii]]
          cols <- seq(mask)
          exY <- summary[traces2average[[iii]], 'holding current (pA)'][cols[mask == 1], drop = FALSE]
          exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
          points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
        }
      }

      Y <- list(summary[, 'holding potential (mV)'])
      if (any(is.na(graph_settings$ylim2))) {
        min_y <- floor(min(summary[, 'holding potential (mV)'], na.rm=TRUE) / 10) * 10
        max_y <- ceiling(max(summary[, 'holding potential (mV)'], na.rm=TRUE) / 10) * 10
        ylim3 <- c(min_y-10, max_y+10)
      }
      plot3(X, Y, xlim=xlim, ylim=ylim3, xlab=xlab, ylab='holding potential (mV)', cex.points=cex.points, 
            title='', colors=colors, ds=1, xmajor_tick=xmajor_tick, ymajor_tick=10, 
            show_spline=FALSE)

      if (length(traces2average) > 0) {
        for (iii in seq_along(split_include)) {
          mask <- split_include[[iii]]
          cols <- seq(mask)
          exY <- summary[traces2average[[iii]], 'holding potential (mV)'][cols[mask == 1], drop = FALSE]
          exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
          points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
        }
      }
      layout(1) 

      if (!force_batch) {
        user_input <- readline("Repeat trace selection? [y/n]: ")
        if (tolower(user_input) %in% c('n', 'no')) {
          break
        }
      } else {
        message("Batch mode: close all plot windows to continue…")
        flush.console()
        # wait until the user has closed every graphics device
        while(length(dev.list()) > 0) {
          Sys.sleep(0.1)
        }
        break
      }
    }
  }

  if (save) {
    wd <- getwd()

    if (any(is.na(graph_settings$ylim1))) {
      max_y <- ceiling(max(summary[, 'peak amplitude (pA)'], na.rm=TRUE) / 100) * 100
      ylim1 <- c(0, max_y)
    }

    # Reverse the order of levels
    setwd(svg_path)
    for (lvl in rev(levels)) {     
      x <- single_examples[['time']]
      if (any(is.na(xlim1))) {
        xlim1 <- c(100*floor(min(x)/100), 100*ceiling(max(x)/100))
      }
      lvl2 <- gsub("[^A-Za-z0-9._-]", "_", lvl)
      filename <- paste0(lvl2, '_eg.svg')
      single_egs(x=x, y=single_examples[[lvl]] , ylim=c(-ylim1[2],20), xlim=xlim1, 
        color='darkgrey', width=width, height=height, filename=filename, save=save)
    }
    setwd(wd)

    # main plot
    dev.new(width=width, height=height, noRStudioGD=TRUE)
    layout(matrix(1:3, ncol = 1), heights = c(3, 2, 2))

    par(mar = c(3.5, 6.5, 2, 1), mgp = c(2.5, 0.5, 0))

    Y <- list(summary[, 'peak amplitude (pA)'])
    plot3(X, Y, xlim=xlim, ylim=ylim1, xlab='', ylab=ylab, cex.points=cex.points, 
          title= as.character(abf_name), colors=colors, ds=1, xmajor_tick=xmajor_tick, 
          show_xaxes=FALSE, ymajor_tick=ymajor_tick, show_spline=FALSE)

    if (length(traces2average) > 0) {
      for (iii in seq_along(split_include)) {
        mask <- split_include[[iii]]
        cols <- seq(mask)
        exY <- summary[traces2average[[iii]], 'peak amplitude (pA)'][cols[mask == 1], drop = FALSE]
        exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
        points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
      }
    }

    if (!is.null(drug_application)) {
      for (i in seq_along(drug_application)) {
        drug_box(
          box_times = c(drug_application[i], tmax),
          box_col = if (i %% 2 == 1) 'darkgrey' else 'lightgrey',
          box_border_col = NA,
          rel_box_height = 0.025,
          rel_y = 1 - (i - 1) * 0.1
        )
      }
    }


    Y <- list(summary[, 'holding current (pA)'])
    if (any(is.na(graph_settings$ylim2))) {
      min_y <- floor(min(summary[, 'holding current (pA)'], na.rm=TRUE) / 100) * 100
      ylim2 <- c(min_y, 0)
    }
    plot3(X, Y, xlim=xlim, ylim=ylim2, xlab='', ylab='holding current (pA)', cex.points=cex.points, 
          title='', colors=colors, ds=1, xmajor_tick=xmajor_tick, 
          ymajor_tick=ymajor_tick, show_xaxes=FALSE, show_spline=FALSE)

    if (length(traces2average) > 0) {
      for (iii in seq_along(split_include)) {
        mask <- split_include[[iii]]
        cols <- seq(mask)
        exY <- summary[traces2average[[iii]], 'holding current (pA)'][cols[mask == 1], drop = FALSE]
        exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
        points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
      }
    }

    Y <- list(summary[, 'holding potential (mV)'])
    if (any(is.na(graph_settings$ylim2))) {
      min_y <- floor(min(summary[, 'holding potential (mV)'], na.rm=TRUE) / 10) * 10
      max_y <- ceiling(max(summary[, 'holding potential (mV)'], na.rm=TRUE) / 10) * 10
      ylim3 <- c(min_y-10, max_y+10)
    }
    plot3(X, Y, xlim=xlim, ylim=ylim3, xlab=xlab, ylab='holding potential (mV)', cex.points=cex.points, 
          title='', colors=colors, ds=1, xmajor_tick=xmajor_tick, ymajor_tick=10, 
          show_spline=FALSE)

    if (length(traces2average) > 0) {
      for (iii in seq_along(split_include)) {
        mask <- split_include[[iii]]
        cols <- seq(mask)
        exY <- summary[traces2average[[iii]], 'holding potential (mV)'][cols[mask == 1], drop = FALSE]
        exX <- summary[traces2average[[iii]], 'time (minutes)'][cols[mask == 1], drop = FALSE]
        points(exX, exY, col='darkgrey', pch=16, cex=cex.points)
      }
    }
    save_graph(svg_path=svg_path, filename=paste0(sub('\\.abf$', '', abf_file), '_PSCs_vs_time.svg'), width=width, height=height, bg='transparent')
    layout(1) 
  }

  levs <- rep(seq_along(traces2average), times = sapply(traces2average, length))
  traces_matrix <- cbind(level=levs, trace=do.call(c, traces2average), include=do.call(c, split_include))

  metadata <- c(
    'dt (ms)' = dt,
    'stimulation time' = stimulation_time,
    'baseline' = baseline, 
    'n' = length(traces2average[[1]])      
  )

  metadata_abf <- extract_metadata(abf_dataset)
  metadata_abf$path <- sub(file_path2, '', metadata_abf$path)

  # create trace selection similar to shiny app (but no rejected column)
  trace_selection <- data.frame(Level = character(), Accepted = character(), 
                                stringsAsFactors = FALSE)
  for (i in seq_along(traces2average)) {
    if (length(traces2average[[i]]) > 0) {
      accepted_indices <- traces2average[[i]][split_include[[i]] == 1]
      accepted_str <- if (length(accepted_indices) > 0) {
        paste(accepted_indices, collapse = ", ")
      } else {
        ""
      }
      trace_selection <- rbind(trace_selection, 
        data.frame(Level = levels[i], Accepted = accepted_str))
    }
  }

  data_list <- list(
    'summary' = summary,
    'raw PSC data' = I_data,
    'baseline corrected' = I_data2,  # CHANGE: 'baseline subtracted PSC data' → 'baseline corrected'
    'trace selection' = trace_selection,  # CHANGE: 'traces2average' with new format
    'single examples' = single_examples,
    'metadata' = t(as.matrix(metadata))
  )

  if (save){

    list2excel(data_list, paste0(sub('\\.abf$', '', abf_file), '.xlsx'), wd=xlsx_path)

    wb <- createWorkbook()
    center_style <- createStyle(halign = 'center', valign = 'center')

    general_df <- do.call(rbind, lapply(names(metadata_abf)[1:5], function(nm) {
      val <- metadata_abf[[nm]]
      value_str <- if (length(val) > 1) paste(val, collapse = ', ') else as.character(val)
      data.frame(Key = nm, Value = value_str, stringsAsFactors = FALSE)
    }))
    addWorksheet(wb, 'General')
    writeData(wb, 'General', general_df)
    addStyle(wb, sheet = 'General', style = center_style,
      rows = 1:(nrow(general_df) + 1), cols = 1:2, gridExpand = TRUE)

    for (nm in names(metadata_abf)[-c(1:5)]) {
      element <- metadata_abf[[nm]]
      if (is.list(element)) {
        df <- flatten_one(element, parent_key = nm)
      } else {
        value_str <- if (length(element) > 1) paste(element, collapse = ', ') else as.character(element)
        df <- data.frame(Key = nm, Value = value_str, stringsAsFactors = FALSE)
      }
      if (!is.null(df) && nrow(df) > 0 && ncol(df) > 0) {
        sheet_name <- substr(nm, 1, 31)
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, df)
        addStyle(wb, sheet = sheet_name, style = center_style,
          rows = 1:(nrow(df) + 1), cols = 1:ncol(df), gridExpand = TRUE)
      }
    }
    saveWorkbook(wb, file = file.path(xlsx_path, paste0(sub('\\.abf$', '', abf_file), '_metadata.xlsx')), overwrite = TRUE)
    graphics.off()
  }
  return(data_list)

}


# helper functions 
# time_stamp <- function(abf_out){
#   hdr <- abf_out$header
#   # parse the date
#   d <- as.character(hdr$uFileStartDate)
#   date <- as.Date(d, format="%Y%m%d")

#   # parse the time
#   ms  <- hdr$uFileStartTimeMS
#   secs <- ms / 1000
#   h <- floor(secs / 3600);       secs <- secs %% 3600
#   m <- floor(secs / 60);         s   <- secs %% 60

#   start_time <- as.POSIXct(sprintf("%s %02d:%02d:%06.3f",
#                                    date, h, m, s),
#                            tz = Sys.timezone())
#   start_time
#   print(start_time, usetz = FALSE)
# }

# time_stamp <- function(abf_out, timezone = "America/Chicago"){
#   hdr <- abf_out$header
#   d <- as.character(hdr$uFileStartDate)
#   date <- as.Date(d, format="%Y%m%d")
#   ms  <- hdr$uFileStartTimeMS
#   secs <- ms / 1000
#   h <- floor(secs / 3600);       secs <- secs %% 3600
#   m <- floor(secs / 60);         s   <- secs %% 60
  
#   # Create POSIXct in specified timezone
#   start_time_local <- as.POSIXct(sprintf("%s %02d:%02d:%02.0f", date, h, m, s),
#                                  tz = timezone)
  
#   # Return formatted string (keeps timezone info)
#   return(format(start_time_local, "%Y-%m-%d %H:%M:%S"))
# }

# Axon Clampex names files based on the session start date (when you started recording that day)
# uFileStartTimeMS gives milliseconds since midnight on that session start date
# If the session spans multiple days, uFileStartTimeMS can exceed 86400000 (24 hours worth of ms)
# This version of time_stamp() handles day overflow to return the correct wall-clock date/time

time_stamp <- function(abf_out, timezone = "America/Chicago"){
  hdr <- abf_out$header
  d <- as.character(hdr$uFileStartDate)
  date <- as.Date(d, format="%Y%m%d")
  ms  <- hdr$uFileStartTimeMS
  secs <- ms / 1000
  h <- floor(secs / 3600);       secs <- secs %% 3600
  m <- floor(secs / 60);         s   <- secs %% 60
  
  # Handle day overflow (sessions spanning midnight)
  extra_days <- h %/% 24
  h <- h %% 24
  date <- date + extra_days
  
  # Create POSIXct in specified timezone
  start_time_local <- as.POSIXct(sprintf("%s %02d:%02d:%02.0f", date, h, m, s),
                                 tz = timezone)
  
  # Return formatted string (keeps timezone info)
  return(format(start_time_local, "%Y-%m-%d %H:%M:%S"))
}


combine_abf_headers <- function(headers) {
  fields   <- names(headers[[1]])
  combined <- list()

  concat_sweeps <- function(v1, v2) {
    # round and coerce to integer
    v1_int <- as.integer(round(v1))
    v2_int <- as.integer(round(v2))

    # infer the uniform sweep interval
    interval <- v1_int[2] - v1_int[1]

    # compute the offset (last start of first recording)
    offset <- max(v1_int)

    # shift *all* v2 starts by offset+interval multiplicatively so
    # first is offset+interval, second is offset+2*interval, etc.
    v2_shifted <- offset + seq_along(v2_int) * interval

    # combine lengths: length(v1_int) + length(v2_int)
    c(v1_int, as.integer(v2_shifted))
  }

  for (f in fields) {
    vals <- lapply(headers, `[[`, f)

    # Scalar numeric fields: sum if not all identical
    if (all(sapply(vals, is.numeric)) && all(sapply(vals, length) == 1)) {
      if (all(vapply(vals[-1], identical, logical(1), vals[[1]]))) {
        combined[[f]] <- vals[[1]]
      } else {
        combined[[f]] <- sum(unlist(vals))
      }

    # sweepStartInPts: special continuous concatenation
    } else if (f == "sweepStartInPts") {
      v1 <- vals[[1]]
      v2 <- vals[[2]]
      combined[[f]] <- concat_sweeps(v1, v2)

    # Other numeric vectors (e.g. telegraph gains): keep identical vector
    } else if (all(sapply(vals, is.numeric)) &&
               all(vapply(vals, length, integer(1)) > 1) &&
               all(vapply(vals[-1], identical, logical(1), vals[[1]]))) {
      combined[[f]] <- vals[[1]]

    # 4) Fallback: store the list of all values
    } else {
      combined[[f]] <- vals
    }
  }

  combined
}

readABFs <- function(abf_files){

  if (length(abf_files)==1){
    out <- readABF(abf_files[1])
  }else{
    out_list <- lapply(abf_files , readABF)
    names(out_list) <- abf_files
    
    # create new abf style output
    master_data <- unlist(lapply(out_list, `[[`, "data"),recursive = FALSE)

    # confirm out_list entries share the same basic metadata
    same_samplingIntervalInSec <- all(sapply(out_list, function(x) 
      identical(x$samplingIntervalInSec, out_list[[1]]$samplingIntervalInSec)
    ))
    same_channelNames <- all(sapply(out_list, function(x) 
      identical(x$channelNames, out_list[[1]]$channelNames)
    ))
    same_channelUnits <- all(sapply(out_list, function(x) 
      identical(x$channelUnits, out_list[[1]]$channelUnits)
    ))
    same_nTelegraphEnable <- all(sapply(out_list, function(x) 
      identical(x$header$nTelegraphEnable, out_list[[1]]$header$nTelegraphEnable)
    ))
    same_fTelegraphAdditGain <- all(sapply(out_list, function(x) 
      identical(x$header$fTelegraphAdditGain, out_list[[1]]$header$fTelegraphAdditGain)
    ))
    same_fInstrumentScaleFactor <- all(sapply(out_list, function(x) 
      identical(x$header$fInstrumentScaleFactor, out_list[[1]]$header$fInstrumentScaleFactor)
    ))

    # If all checks pass, build the master object
    if (all(
      same_channelNames,
      same_channelUnits,
      same_nTelegraphEnable,
      same_fTelegraphAdditGain,
      same_fInstrumentScaleFactor
    )) {
      out <- out_list[[1]]
      out$data <- master_data
      headers    <- lapply(out_list, `[[`, "header")
      out$header <- combine_abf_headers(headers)
    } else {
      stop("ABF metadata not uniform across all files")
    }
  }
  out
}

focus_console <- function() {
  script <- paste0(
    'tell application "System Events"\n',
    '  set frontmost of process "R" to true\n',
    'end tell'
  )
  system2("osascript", c("-e", shQuote(script)))
}

select_traces <- function(all_traces, dt,
                          ymajor_tick = 100, xmajor_tick = 200,
                          amp = NULL, hc = NULL,
                          mgp = c(2, 0.5, 0),
                          height = 5, width = 5, cex = 0.7,
                          default_include = 1,
                          force_batch = FALSE) {
  time <- seq(0, nrow(all_traces) - 1) * dt
  n <- ncol(all_traces)
  include <- integer(n)

  for (ii in seq_len(n)) {
    y <- all_traces[, ii]
    y_range <- range(y, na.rm = TRUE)
    ylim <- c(floor(y_range[1]/100)*100, ceiling(y_range[2]/100)*100)
    xlim <- c(0, ceiling(max(time)/100)*100)

    if (force_batch || !interactive()) {
      # Batch mode: no plotting, auto-include
      include[ii] <- default_include
      next
    }

    # Interactive plotting
    main <- if (is.null(amp)) {
      paste('Trace', ii)
    } else {
      paste0('Trace', ii, '; amp=', round(amp[ii], 2),
             ' pA; hc=', round(hc[ii], 2), ' pA')
    }

    par(mgp = mgp)
    plot(time, y, type = 'l', main = main,
         ylim = ylim, xlim = xlim,
         xlab = 'time (ms)', ylab = 'pA', axes = FALSE,
         cex.main = 10/7*cex, cex.lab = cex, cex.axis = cex,
         xaxs = 'i', yaxs = 'i', bty = 'n', col = '#A9A9A9')

    axis(1, at = seq(xlim[1], xlim[2], by = xmajor_tick),
         tcl = -0.1, cex.axis = cex)
    axis(2, at = seq(ylim[1], ylim[2], by = ymajor_tick),
         tcl = -0.1, cex.axis = cex, las = 1)

    # Prompt on the same line
    flush.console()
    response <- readline(
      prompt = paste0("include trace ", ii, " ? (1 = Yes, 0 = No): ")
    )
    while (!(response %in% c("0", "1"))) {
      flush.console()
      response <- readline(
        prompt = "please enter 1 to include, 0 to exclude: "
      )
    }
    include[ii] <- as.integer(response)
  }

  include
}

extract_metadata <- function(abf_dataset) {
    list(
      path                  = abf_dataset$path,
      formatVersion         = abf_dataset$formatVersion,
      channelNames          = abf_dataset$channelNames,
      channelUnits          = abf_dataset$channelUnits,
      samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
      header                = abf_dataset$header,
      tags                  = abf_dataset$tags,
      sections              = abf_dataset$sections
    )
  }

extract_response <- function(y, dt, stimulation_time, baseline){
  
  idx1 <- (stimulation_time - baseline) / dt
  idx2 <- baseline / dt

  y1 <- y[idx1:length(y)]
  y1 <- y1 - mean(y1[0:idx2])

  y1
}

# if ms and pA then output would be fC so 1e3 corrects to pC
trap_fun <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / 1e3
}

peak.fun2 <- function(y, dt, stimulation_time, baseline, detection_window=200, smooth=5){
  
  idx1 <- (stimulation_time - baseline) / dt
  idx2 <- baseline / dt
  idx3 <-  detection_window / dt

  y1 <- y[idx1:length(y)]
  bl <- mean(y1[0:idx2])

  y1 <- y1 - bl
  peak <- moving_avg(y1[idx2:(idx2 + idx3)], n = smooth)

  y2 <- y1[idx2:length(y1)]
  x <- (1:length(y2))*dt
  area <- trap_fun(x, y2)

  return(list(response=y1, peak=peak, charge=area, baseline=bl))
}

# helper to flatten abf metadata before export
flatten_one <- function(x, parent_key = '') {
  if (is.atomic(x) || is.null(x)) {
    key <- parent_key
    val <- if (length(x) > 1) paste(x, collapse = ', ') else as.character(x)
    return(data.frame(Key = key, Value = val, stringsAsFactors = FALSE))
  }

  if (is.list(x)) {
    do.call(rbind, Map(function(nm, val) {
      flatten_one(val, paste0(parent_key, if (parent_key != '') '$', nm))
    }, names(x), x))
  } else {
    return(data.frame(Key = parent_key, Value = as.character(x), stringsAsFactors = FALSE))
  }
}

moving_avg <- function(y, n = 5) {
  sign <- sign_fun(y)
  y <- y * sign
  y_length <- length(y)
  result <- rep(NA, y_length)
  
  for (i in 1:y_length) {
    # Determine the start and end indices for the window
    start_idx <- max(1, i - floor(n / 2))
    end_idx <- min(y_length, i + floor(n / 2))
    
    # Calculate the mean for the current window
    result[i] <- mean(y[start_idx:end_idx], na.rm = TRUE)
  }
  
  return(max(result) * sign)
}

sign_fun <- function(y, direction_method = c('smooth', 'regression', 'cumsum'), k = 5) {
  
  method <- match.arg(direction_method)
  
  # Check if the length of y is sufficient for processing
  if (length(y) < 10) {
    stop("The length of y must be at least 10.")
  }
  
  n <- length(y)
  peak_value <- NA
  
  if (method == 'smooth') {
    # Calculate the smoothed signal using a simple moving average
    smoothed_signal <- rep(NA, n)
    for (i in 1:(n - k + 1)) {
      smoothed_signal[i + floor(k/2)] <- mean(y[i:(i + k - 1)])
    }
    # Remove NA values before finding the indices of the maximum and minimum values
    valid_indices <- which(!is.na(smoothed_signal))
    max_idx <- valid_indices[which.max(smoothed_signal[valid_indices])]
    min_idx <- valid_indices[which.min(smoothed_signal[valid_indices])]
    
    # Determine the peak value based on the magnitude comparison
    if (abs(smoothed_signal[max_idx]) > abs(smoothed_signal[min_idx])) {
      peak_value <- smoothed_signal[max_idx]
    } else {
      peak_value <- smoothed_signal[min_idx]
    }
    
  } else if (method == 'diff') {
    # Calculate the differences of the signal
    diff_signal <- diff(y)
    # Identify both the maximum and minimum differences
    max_diff <- max(diff_signal, na.rm = TRUE)
    min_diff <- min(diff_signal, na.rm = TRUE)
    
    # Determine the peak value considering the direction
    peak_value <- ifelse(abs(max_diff) > abs(min_diff), max_diff, min_diff)
     
  } else if (method == 'regression') {
    # Fit a quadratic regression to the signal
    x <- 1:n
    model <- lm(y ~ I(x^2) + x)
    
    # Use the fitted values to find the peak
    fitted_values <- predict(model)
    max_idx <- which.max(fitted_values)
    min_idx <- which.min(fitted_values)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(fitted_values[max_idx]) > abs(fitted_values[min_idx])) {
      peak_value <- fitted_values[max_idx]
    } else {
      peak_value <- fitted_values[min_idx]
    }
    
  } else if (method == 'cumsum') {
    # Calculate the cumulative sum of the signal
    cumsum_signal <- cumsum(y)
    # Find the indices of the maximum and minimum values
    max_idx <- which.max(cumsum_signal)
    min_idx <- which.min(cumsum_signal)
    
    # Determine the peak value based on the magnitude comparison
    if (abs(cumsum_signal[max_idx]) > abs(cumsum_signal[min_idx])) {
      peak_value <- cumsum_signal[max_idx]
    } else {
      peak_value <- cumsum_signal[min_idx]
    }
  }
  
  # Determine if the peak is positive or negative
  peak_direction <- ifelse(peak_value > 0, 1, -1)
  
  # Return the sign of the peak direction
  return(peak_direction)
}

plot2 <- function(x_list, y_list, colors=NULL, xlim=c(0, 275), ylim=c(0, 70), lwd=1, 
  xlab='distance from soma center (\u00B5m)', ylab='PSP amplitude (mV)', height=3, width=3, 
  bg='transparent', h=NULL, xmajor_tick=50, ymajor_tick=20, ds=3, title=NULL, cex=0.7, 
  cex.points=0.3, mgp=c(2, 0.5, 0), xaxes_offset=0.025, yaxes_offset=0.02, spar=NULL, 
  show_spline=TRUE, filename='plot2.svg', save=FALSE, logx=FALSE) {

  # Remove zero or negative x values if logx is TRUE
  if (logx) {
    valid_data <- mapply(function(x, y) {
      valid_idx <- x > 0
      list(x = x[valid_idx], y = y[valid_idx])
    }, x_list, y_list, SIMPLIFY = FALSE)
    
    x_list <- lapply(valid_data, `[[`, "x")
    y_list <- lapply(valid_data, `[[`, "y")
  }
  
  # Extend axis range
  x_range <- diff(xlim)
  if (!logx) {
    xlim <- c(xlim[1] - xaxes_offset * x_range, xlim[2] + xaxes_offset * x_range)
  } else {
    xlim <- c(max(xlim[1] * (1 - xaxes_offset), .1), xlim[2] * (1 + xaxes_offset))
  }
  
  y_range <- diff(ylim)
  ylim <- c(ylim[1] - yaxes_offset * y_range, ylim[2] + yaxes_offset * y_range)
  
  # Open SVG device
  n <- length(y_list)
  if (is.null(colors)) colors <- hex_palette2(n=n, reverse=FALSE)
  
  if (save) {
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  # Set graphical parameters
  par(mgp=mgp)

  # Create plot
  plot(NA, type='n', xlim=xlim, ylim=ylim, bty='n', main=title, cex.main=cex, axes=FALSE, 
       cex.lab=cex, cex.axis=cex, xaxs='i', yaxs='i', xlab=xlab, ylab=ylab, log=ifelse(logx, "x", ""))
  
  # Loop through the y_list and plot each smooth spline
  for (ii in seq_along(y_list)) {
    x_points <- x_list[[ii]]
    y_points <- y_list[[ii]]
    n_points <- length(x_points)
    
    # Downsample
    if (n_points > 2) {
      indices <- c(1, seq(2, n_points - 1, by=ds), n_points)
    } else {
      indices <- 1:n_points
    }
    
    # Plot the points
    points(x_points[indices], y_points[indices], col=colors[ii], pch=16, cex=cex.points)
    
    # Smooth spline
    if (show_spline){
      if (length(x_points) > 1) {
        if (logx) {
          # Transform x_points to log10 space for spline fitting
          log_x_points <- log10(x_points)
          spline_fit <- smooth.spline(log_x_points, y_points, spar=spar)
          xspline <- seq(min(log_x_points), max(log_x_points), length.out=1000)
          yspline <- predict(spline_fit, xspline)$y
          xspline <- 10^xspline  # Transform back to linear space for plotting
          # print(xspline[which.min(yspline)])
        } else {
          spline_fit <- smooth.spline(x_points, y_points, spar=spar)
          xspline <- seq(min(x_points), max(x_points), length.out=1000)
          yspline <- predict(spline_fit, xspline)$y
          # print(xspline[which.min(yspline)])
        }
        # Plot the smooth spline
        lines(xspline, yspline, col=colors[ii], lwd=lwd, lty=3)
      }
    }
  }
  
  # Add vertical lines if h is specified
  if (!is.null(h)) {
    abline(v = h, col = 'darkgray', lty=3, lwd=lwd)
  }

  # Add x-axis ticks (logarithmic when logx = TRUE)
  if (!logx) {
    axis(1, at=seq(xlim[1] + xaxes_offset * x_range, xlim[2] - xaxes_offset * x_range, by=xmajor_tick), tcl=-0.2, cex.axis=cex) 
    axis(1, at=seq(xlim[1] + xaxes_offset * x_range, xlim[2] - xaxes_offset * x_range, by=xmajor_tick / 2), labels=FALSE, tcl=-0.1) 
  } else {
    # Logarithmic major ticks
    log_range <- log10(xlim)  # Log10 of the xlim range

    # Major ticks: 0.1, 1, 3, 10, 30...
    major_ticks <- 10^seq(floor(log_range[1]), ceiling(log_range[2]), by=1)  # Decade ticks (1, 10, 100...)
    major_ticks <- major_ticks[major_ticks >= xlim[1] & major_ticks <= xlim[2]]  # Filter within xlim

    # Minor ticks: 2, 5 between major ticks
    minor_ticks <- seq(1:9) %o% 10^(floor(log_range[1]):ceiling(log_range[2]))  # Minor ticks at each decade
    minor_ticks <- as.vector(minor_ticks)
    minor_ticks <- minor_ticks[minor_ticks >= xlim[1] & minor_ticks <= xlim[2]]  # Filter within xlim

    # Add major ticks and labels
    axis(1, at=major_ticks, labels=major_ticks, tcl=-0.2, cex.axis=cex)

    # Add minor ticks (no labels)
    axis(1, at=minor_ticks, labels=NA, tcl=-0.1)

  }
  
  # Add y-axis ticks
  axis(2, at=seq(ylim[1] + yaxes_offset * y_range, ylim[2] - yaxes_offset * y_range, by=ymajor_tick), las=2, tcl=-0.2, cex.axis=cex) 
  axis(2, at=seq(ylim[1] + yaxes_offset * y_range, ylim[2] - yaxes_offset * y_range, by=ymajor_tick / 2), labels=FALSE, tcl=-0.1) 
  
  if (save) dev.off()
}

plot3 <- function(x_list, y_list, colors=NULL, xlim=c(0, 275), ylim=c(0, 70), lwd=1, 
  xlab='distance from soma center (\u00B5m)', ylab='PSP amplitude (mV)', 
  bg='transparent', h=NULL, xmajor_tick=50, ymajor_tick=20, ds=3, title=NULL, cex=1, 
  cex.points=0.3, xaxes_offset=0.025, yaxes_offset=0.02, spar=NULL, 
  show_xaxes=TRUE, show_spline=TRUE, logx=FALSE) {

  if (logx) {
    valid_data <- mapply(function(x, y) {
      valid_idx <- x > 0
      list(x = x[valid_idx], y = y[valid_idx])
    }, x_list, y_list, SIMPLIFY = FALSE)
    
    x_list <- lapply(valid_data, `[[`, "x")
    y_list <- lapply(valid_data, `[[`, "y")
  }
  
  x_range <- diff(xlim)
  if (!logx) {
    xlim <- c(xlim[1] - xaxes_offset * x_range, xlim[2] + xaxes_offset * x_range)
  } else {
    xlim <- c(max(xlim[1] * (1 - xaxes_offset), .1), xlim[2] * (1 + xaxes_offset))
  }
  
  y_range <- diff(ylim)
  ylim <- c(ylim[1] - yaxes_offset * y_range, ylim[2] + yaxes_offset * y_range)
  
  n <- length(y_list)
  if (is.null(colors)) colors <- hex_palette2(n=n, reverse=FALSE)


  plot(NA, type='n', xlim=xlim, ylim=ylim, bty='n', main=title, cex.main=cex, axes=FALSE, 
       cex.lab=cex, cex.axis=cex, xaxs='i', yaxs='i', xlab=xlab, ylab=ylab, log=ifelse(logx, "x", ""))

  for (ii in seq_along(y_list)) {
    x_points <- x_list[[ii]]
    y_points <- y_list[[ii]]
    n_points <- length(x_points)
    
    if (n_points > 2) {
      indices <- c(1, seq(2, n_points - 1, by=ds), n_points)
    } else {
      indices <- 1:n_points
    }
    
    points(x_points[indices], y_points[indices], col=colors[ii], pch=16, cex=cex.points)
    
    if (show_spline){
      if (length(x_points) > 1) {
        if (logx) {
          log_x_points <- log10(x_points)
          spline_fit <- smooth.spline(log_x_points, y_points, spar=spar)
          xspline <- seq(min(log_x_points), max(log_x_points), length.out=1000)
          yspline <- predict(spline_fit, xspline)$y
          xspline <- 10^xspline
        } else {
          spline_fit <- smooth.spline(x_points, y_points, spar=spar)
          xspline <- seq(min(x_points), max(x_points), length.out=1000)
          yspline <- predict(spline_fit, xspline)$y
        }
        lines(xspline, yspline, col=colors[ii], lwd=lwd, lty=3)
      }
    }
  }

  if (!is.null(h)) {
    abline(v = h, col = 'darkgray', lty=3, lwd=lwd)
  }

  if (show_xaxes){
    if (!logx) {
      axis(1, at=seq(xlim[1] + xaxes_offset * x_range, xlim[2] - xaxes_offset * x_range, by=xmajor_tick), tcl=-0.2, cex.axis=cex) 
      axis(1, at=seq(xlim[1] + xaxes_offset * x_range, xlim[2] - xaxes_offset * x_range, by=xmajor_tick / 2), labels=FALSE, tcl=-0.1) 
    } else {
      log_range <- log10(xlim)
      major_ticks <- 10^seq(floor(log_range[1]), ceiling(log_range[2]), by=1)
      major_ticks <- major_ticks[major_ticks >= xlim[1] & major_ticks <= xlim[2]]
      minor_ticks <- seq(1:9) %o% 10^(floor(log_range[1]):ceiling(log_range[2]))
      minor_ticks <- as.vector(minor_ticks)
      minor_ticks <- minor_ticks[minor_ticks >= xlim[1] & minor_ticks <= xlim[2]]
      axis(1, at=major_ticks, labels=major_ticks, tcl=-0.2, cex.axis=cex)
      axis(1, at=minor_ticks, labels=NA, tcl=-0.1)
    }
  }
  
  axis(2, at=seq(ylim[1] + yaxes_offset * y_range, ylim[2] - yaxes_offset * y_range, by=ymajor_tick), las=2, tcl=-0.2, cex.axis=cex) 
  axis(2, at=seq(ylim[1] + yaxes_offset * y_range, ylim[2] - yaxes_offset * y_range, by=ymajor_tick / 2), labels=FALSE, tcl=-0.1) 
}

# drug application box
drug_box <- function(box_times, box_col='darkgrey', box_border_col=NA, rel_box_height=0.05, rel_y=1){
  usr <- par('usr')               # c(xleft, xright, ybottom, ytop)
  top <- usr[4]                   # top of y range
  top1 <- rel_y*usr[4]            # top of drug box
  heightFrac <- rel_box_height    # box height as fraction of total y range
  yBoxBottom <- top1 - (top - usr[3]) * heightFrac

  rect(xleft = box_times[1], 
       ybottom = yBoxBottom,
       xright = box_times[2],
       ytop = top1,
       col = box_col,  # or your choice of fill color
       border = box_border_col)  # NA or some color as required
}

single_egs <- function(x, y, sign=-1, xlim=NULL, ylim=NULL, lwd=1, show_text=FALSE, height=4, width=2.5, xbar=100, ybar=50,  color='#4C77BB', filename='trace1.svg', bg='transparent', save=FALSE) {
  
  if (save) {
    svg(file=filename, width=width, height=height, bg=bg)
  } else {
    dev.new(width=width, height=height, noRStudioGD=TRUE)
  }

  if (is.null(ylim)) ylim <- if (sign == 1) c(0, max(y)) else c(-max(-y), 0)

  if (is.null(xlim)) xlim <- c(min(x), max(x))
  idx1 <- which.min(abs(x - xlim[1]))
  idx2 <- which.min(abs(x - xlim[2]))

  plot(x[idx1:idx2], y[idx1:idx2], type='l', col=color, xlim=xlim, ylim=ylim, bty='n', lwd=lwd, lty=1, axes=FALSE, frame=FALSE, xlab='', ylab='')

  #  scale bar lengths and ybar position
  ybar_start <- min(ylim) + (max(ylim) - min(ylim)) / 20
  
  # Add scale bars at the bottom right
  x_start <- max(xlim) - xbar - 50
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar
  
  # Draw the scale bars
  segments(x_start, y_start, x_end, y_start, lwd=lwd, col='black')
  segments(x_start, y_start, x_start, y_end, lwd=lwd, col='black')
  
  # Add labels to the scale bars
  if (show_text) {
    text(x = (x_start + x_end) / 2, y = y_start - ybar / 20, labels = paste(xbar, 'ms'), adj = c(0.5, 1))
    text(x = x_start - xbar / 4, y = (y_start + y_end) / 2,  labels = paste(ybar, 'pA'), adj = c(0.5, 0.5), srt = 90)
  }
    
  if (save) {
    dev.off()
  }
}

# can use this function to save any open graph 
save_graph <- function(svg_path, filename='graph1.svg', width=6, height=4, bg="transparent") {
  old_wd <- getwd()
  setwd(svg_path)
  dev.copy(svg, file=filename, width=width, height=height, bg=bg)
  dev.off()
  setwd(old_wd)
}

# save list to excel spreadsheet
list2excel <- function(data_list, file_name, wd = getwd(), center_align = TRUE) {
  # Load the openxlsx library
  library(openxlsx)
  
  # Create a new workbook
  workbook <- createWorkbook()
  
  # Define styles: bold and centered for headers, and centered for data
  header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
  center_style <- createStyle(halign = "center", valign = "center")
  
  # Loop over each element in the list
  for (i in seq_along(data_list)) {
    # Use the name of the list element as the sheet name
    sheet_name <- names(data_list)[i]
    
    # Default to "Sheet1", "Sheet2", etc., if the name is missing
    if (is.null(sheet_name) || sheet_name == "") {
      sheet_name <- paste0("Sheet", i)
    }
    
    # Add a new sheet with the specified name to the workbook
    addWorksheet(workbook, sheet_name)
    
    # Write the data to the sheet
    writeData(workbook, sheet_name, data_list[[i]])
    
    # Apply header style to make headers bold and centered
    addStyle(workbook, sheet_name, style = header_style, rows = 1, cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    
    # Apply center alignment to all cells if center_align is TRUE
    if (center_align) {
      addStyle(workbook, sheet_name, style = center_style, rows = 1:(nrow(data_list[[i]]) + 1), 
               cols = 1:ncol(data_list[[i]]), gridExpand = TRUE)
    }
  }
  
  # Create the full file path
  file_path <- file.path(wd, file_name)
  
  # Save the workbook
  saveWorkbook(workbook, file_path, overwrite = TRUE)
}


amplifier_gain <- function(dataset = NULL, headstage_gain = 0.5, additional_gain = NULL,
                           AD_range = c(-10, 10), AD_bits = 16,
                           dp = 3, tol = 1e-3, VClamp = TRUE) {
  
  if (!is.null(dataset) && !is.null(additional_gain)) {
    amplifier_gain3(dataset = dataset,
                    headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    tol = tol,
                    dp = dp,
                    VClamp = VClamp)
    
  } else if (!is.null(dataset)) {
    amplifier_gain1(dataset = dataset,
                    headstage_gain = headstage_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    dp = dp,
                    tol = tol,
                    VClamp = VClamp)
    
  } else {
    amplifier_gain2(headstage_gain = headstage_gain,
                    additional_gain = additional_gain,
                    AD_range = AD_range,
                    AD_bits = AD_bits,
                    VClamp = VClamp)
  }
}

dpA_fun <- function(dataset, tol=1e-2){
  sapply(1:dim(dataset)[2], function(ii){
    vec <- diff(sort(unique(dataset[,ii])))
    vec <- vec[vec>tol]
    min(vec)
    }
  )
}

amplifier_gain1 <- function(dataset, headstage_gain=0.5, AD_range=c(-10, 10), AD_bits=16, dp=3, tol=1e-3, VClamp=TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  
  if (VClamp) {
    dpA <- dpA_fun(dataset=dataset, tol=tol)
    recording_range <- dpA * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain mV/pA' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dpA,
      check.names = FALSE
    )
    
  } else {
    dV <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    recording_range <- dV * 2^AD_bits
    min_recording <- -recording_range / 2
    max_recording <-  recording_range / 2
    final_gain <- digitiser_range * 1e3 / recording_range
    additional_gain <- final_gain / headstage_gain
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = round(additional_gain, dp),
      'final gain V/V' = round(final_gain, dp),
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dV,
      check.names = FALSE
    )
  }
  
  return(output)
}

amplifier_gain2 <- function(headstage_gain=0.5, additional_gain=20, AD_range=c(-10, 10), AD_bits=16, VClamp=TRUE) {
  
  if (length(additional_gain) == 1 && length(headstage_gain) > 1) {
    additional_gain <- rep(additional_gain, length(headstage_gain))
  }
  
  if (length(headstage_gain) != length(additional_gain)) {
    stop("if 'additional_gain' is a vector, it must be the same length as 'headstage_gain'")
  }
  
  min_A_D <- AD_range[1]
  max_A_D <- AD_range[2]
  digitiser_range <- abs(diff(AD_range))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits
  
  if (VClamp) {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain mV/pA' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain mV/pA' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      check.names = FALSE
    )
  } else {
    output <- data.frame(
      'R GOhms' = headstage_gain,
      'gain V/V' = headstage_gain,
      'additional gain' = additional_gain,
      'final gain V/V' = final_gain,
      'min A-D board V' = rep(min_A_D, length(headstage_gain)),
      'max A-D board V' = rep(max_A_D, length(headstage_gain)),
      'A-D board range V' = rep(digitiser_range, length(headstage_gain)),
      'A-D bits' = rep(AD_bits, length(headstage_gain)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      check.names = FALSE
    )
  }
  return(output)
}

amplifier_gain3 <- function(dataset, headstage_gain = 0.5, additional_gain = 20,
                            AD_range = c(-10, 10), AD_bits = 16,
                            tol = 1e-3, dp=3, VClamp = TRUE) {
  
  digitiser_range <- abs(diff(AD_range))
  min_A_D <- rep(AD_range[1], ncol(dataset))
  max_A_D <- rep(AD_range[2], ncol(dataset))
  final_gain <- headstage_gain * additional_gain
  recording_range <- digitiser_range * 1e3 / final_gain
  min_recording <- -recording_range / 2
  max_recording <-  recording_range / 2
  dUnit <- recording_range / 2^AD_bits  # actual theoretical digitisation
  
  if (VClamp) {
    dpA_actual <- dpA_fun(dataset = dataset, tol = tol)
    n <- dUnit / dpA_actual 
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain mV/pA' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain mV/pA' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording pA' = min_recording,
      'max recording pA' = max_recording,
      'recording range pA' = recording_range,
      'digitisation pA/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
    
  } else {
    dV_actual <- sapply(1:dim(dataset)[2], function(ii) min(diff(sort(unique(dataset[, ii])))) )
    n <- dUnit / dV_actual
    
    output <- data.frame(
      'R GOhms' = rep(headstage_gain, ncol(dataset)),
      'gain V/V' = rep(headstage_gain, ncol(dataset)),
      'additional gain' = rep(additional_gain, ncol(dataset)),
      'final gain V/V' = final_gain,
      'min A-D board V' = min_A_D,
      'max A-D board V' = max_A_D,
      'A-D board range V' = rep(digitiser_range, ncol(dataset)),
      'A-D bits' = rep(AD_bits, ncol(dataset)),
      'min recording mV' = min_recording,
      'max recording mV' = max_recording,
      'recording range mV' = recording_range,
      'digitisation mV/bit' = dUnit,
      'n' = round(n, dp),
      check.names = FALSE
    )
  }
  
  return(output)
}

format_experimental_dict <- function(dict) {
  for (group in names(dict)) {
    
    for (folder in names(dict[[group]])) {
      entry <- dict[[group]][[folder]]
      line <- paste0("    '", folder, "' = list(")

      # traces2average first
      if (!is.null(entry$traces2average)) {
        traces_str <- paste0(
          "traces2average = list(",
          paste(
            sapply(entry$traces2average, function(vec) {
              paste0("c(", paste(vec, collapse = ","), ")")
            }),
            collapse = ", "
          ),
          "),    "
        )
        line <- paste0(line, traces_str)
      }

      # then indices and levels
      line <- paste0(line, "indices = c(", paste(entry$indices, collapse = ","), "),    ")
      line <- paste0(line, "levels = c('", paste(entry$levels, collapse = "','"), "')")

      line <- paste0(line, "),")
      cat(line, "\n")
    }
    cat("  ),\n")
  }
}


load_data2 <- function(wd, name, header = TRUE) {
  # Create the file path
  file_path <- file.path(wd, paste0(name, '.', 'xlsx'))
  
  # Load the Excel file
  workbook <- openxlsx::loadWorkbook(file_path)
  
  # Get the sheet names
  sheet_names <- openxlsx::getSheetNames(file_path)
  
  # Initialize an empty list to store each sheet's data
  data_list <- list()
  
  # Loop through each sheet and read the data into the list
  for (sheet in sheet_names) {
    data_list[[sheet]] <- openxlsx::read.xlsx(file_path, sheet = sheet, colNames = header)
  }
  
  return(data_list)
}

digitisation_fun <- function(x) {
  diffs <- abs(diff(x))
  nz_diffs <- diffs[diffs > 0]
  if (length(nz_diffs) == 0) return(NA_real_)
  min(nz_diffs)
}

get_all_digitisation <- function(x_list, col) {
  sapply(x_list, function(mat) digitisation_fun(mat[, col]))
}


get_adc_scaling_summary <- function(meta, AD_range = c(-10, 10), AD_bits = 16) {
  adc_info <- meta$sections$ADCsec
  units <- meta$channelUnits
  digitisation <- meta$digitisation
  n <- min(2, length(adc_info))

  result <- lapply(seq_len(n), function(i) {
    ch <- adc_info[[i]]
    unit <- units[i]
    scale <- ch$fInstrumentScaleFactor * 1000
    gain <- ch$fTelegraphAdditGain
    total_gain <- gain * scale

    V_range <- diff(AD_range)
    LSB <- V_range / (2^AD_bits)
    AD_max <- max(AD_range) / (total_gain / 1000)
    AD_min <- min(AD_range) / (total_gain / 1000)
    AD_LSB <- LSB / (total_gain / 1000)

    # Extract trace digitisation values
    trace_vals <- digitisation[[paste0('col', i)]]
    trace_vals <- trace_vals[!is.na(trace_vals)]
    trace_digitisation <- if (length(trace_vals) > 0 && all(abs(trace_vals - trace_vals[1]) < .Machine$double.eps)) {
      trace_vals[1]
    } else {
      NA_real_
    }

    df <- data.frame(
      channel = ch$nADCNum,
      units = unit,
      'telegraph enabled' = as.logical(ch$nTelegraphEnable),
      'scale factor (mV/unit)' = scale,
      'additional gain' = gain,
      'total gain (mV/unit)' = total_gain,
      'gain units' = paste0('mV/', unit),
      'AD min' = AD_min,
      'AD max' = AD_max,
      'digitisation' = AD_LSB,
      'trace digitisation' = trace_digitisation,
      recorded = (ch$nTelegraphEnable == 1 && gain > 0 && scale > 0),
      check.names = FALSE
    )
    df
  })

  do.call(rbind, result)
}


# FUNCTIONS FOR NWB EXPORT


# Extract metadata from ABF dataset
extract_metadata <- function(abf_dataset) {
  list(
    path                  = abf_dataset$path,
    formatVersion         = abf_dataset$formatVersion,
    channelNames          = abf_dataset$channelNames,
    channelUnits          = abf_dataset$channelUnits,
    samplingIntervalInSec = abf_dataset$samplingIntervalInSec,
    header                = abf_dataset$header,
    tags                  = abf_dataset$tags,
    sections              = abf_dataset$sections
  )
}



# MODULAR VERSION OF analyseABF() WITH NWB EXPORT
# UI HELPER FUNCTIONS
create_dark_mode_css <- function() {
  tags$head(
    tags$style(HTML("
      @media (prefers-color-scheme: dark) {
        body { background-color: #1e1e1e; color: #e0e0e0; }
        .well { background-color: #2d2d2d; border-color: #444; }
        .form-control { background-color: #2d2d2d; color: #c0c0c0 !important; border: 1px solid #555 !important; }
        input[type='number'], input[type='text'] { background-color: #2d2d2d !important; color: #c0c0c0 !important; }
        .selectize-input, .selectize-dropdown { background-color: #ffffff !important; color: #666666 !important; }
        .btn { background-color: #3d3d3d; color: #ffffff; border-color: #555; }
        .btn-primary, .action-button { background-color: #3c8dbc; color: #ffffff; }
        pre, code { background-color: #1a1a1a; color: #f0f0f0; }
      }
    "))
  )
}

create_main_settings_ui <- function() {
  tabPanel("Main Settings",
    numericInput('baseline', 'Baseline (ms):', 100, min = 0),
    numericInput('stimulation', 'Stimulation Time (ms):', value = 150),
    checkboxInput('autoDetectStim', 'Auto-detect from TTL', TRUE),
    hr(),
    textInput('levels', 'Levels (comma-separated):', 'control,drug'),
    textInput('drugApplication', 'Drug Application Times (comma-separated):', ''),
    hr(),
    numericInput('pscChannel', 'PSC Channel:', 1, min = 1, max = 10),
    numericInput('hpChannel', 'Holding Potential Channel:', 2, min = 1, max = 10),
    numericInput('ttlChannel', 'TTL Channel (optional):', 3, min = 1, max = 10),
    hr(),
    verbatimTextOutput('fileInfo')
  )
}

create_trace_selection_ui <- function() {
  tabPanel("Trace Selection",
    helpText("Define traces to average for each level"),
    uiOutput('traceSelectionUI'),
    hr(),
    verbatimTextOutput('tracesInfo')
  )
}

create_graph_settings_ui <- function() {
  tabPanel("Graph Settings",
    numericInput('width', 'Plot Width:', 6, min = 1, max = 20),
    numericInput('height', 'Plot Height:', 8, min = 1, max = 20),
    hr(),
    h5("Peak Amplitude Plot"),
    numericInput('xmajor_tick_amp', 'X Tick:', 5, min = 1),
    numericInput('ymajor_tick_amp', 'Y Tick:', 100, min = 1),
    hr(),
    h5("Holding Current Plot"),
    numericInput('ymajor_tick_hc', 'Y Tick:', 100, min = 1),
    hr(),
    h5("Holding Potential Plot"),
    numericInput('xmajor_tick_hp', 'X Tick:', 5, min = 1),
    numericInput('ymajor_tick_hp', 'Y Tick:', 10, min = 1),
    hr(),
    h5("Appearance"),
    numericInput('cex_points', 'Point Size:', 1.5, min = 0.1, max = 3, step = 0.1),
    numericInput('lwd_graph', 'Line Thickness:', 1.5, min = 0.5, max = 5, step = 0.1),
    numericInput('cex_labels', 'Label Size:', 1.2, min = 0.5, max = 3, step = 0.1),
    numericInput('cex_axis', 'Axis Text Size:', 1.6, min = 0.5, max = 3, step = 0.1),
    hr(),
    textInput('xlab', 'X Label:', 'time (minutes)'),
    textInput('ylab', 'Y Label:', '|PSC| (pA)'),
    textInput('color', 'Plot Color:', '#CD5C5C')
  )
}

create_nwb_settings_ui <- function() {
  tabPanel("NWB Export",
    h5("Subject Information"),
    textInput('nwb_subject_id', 'Subject ID:', 'Mouse001'),
    textInput('nwb_species', 'Species:', 'Mus musculus'),
    textInput('nwb_age', 'Age:', 'P30D'),
    selectInput('nwb_sex', 'Sex:', choices = c('M', 'F', 'U'), selected = 'M'),
    textInput('nwb_genotype', 'Genotype:', 'WT'),
    textInput('nwb_cross', 'Cross:', ''),
    textInput('nwb_virus', 'Virus:', ''),
    textInput('nwb_virus_injection_site', 'Virus Injection Site:', ''),
    hr(),
    h5("Session Information"),
    textInput('nwb_session_id', 'Session ID:', 's01'),
    textInput('nwb_task', 'Task:', 'whc'),
    hr(),
    h5("Experimental Details"),    
    textInput('nwb_location', 'Recording Location:', 'Striatum'),
    textInput('nwb_electrode_resistance', 'Electrode Resistance:', '3-5 MOhm'),
    hr(),
    h5("Lab Information"),
    textInput('nwb_experimenter', 'Experimenter (LastName, FirstName):', 'LastName, FirstName'), 
    textInput('nwb_institution', 'Institution:', 'YourInstitution'),
    textInput('nwb_lab', 'Lab:', 'YourLab'),
    hr(),
    h5("Experiment Description"),                                                          
    tags$textarea(                                                                         
      id = "nwb_experiment_description",                                                  
      class = "form-control",                                                             
      rows = "3",                                                                          
      placeholder = "Brief description of the overall experiment..."                     
    ),                                                                                     
    textInput('nwb_keywords', 'Keywords (comma-separated):', 'striatum,voltage clamp'),   
    hr(),                                                                                  
    h5("Device"),
    textInput('nwb_device_name', 'Device:', 'Multiclamp 700B'),
    hr(),
    tags$textarea(                                                                        
      id = "nwb_device_description",                                                      
      class = "form-control",                                                             
      rows = "2",                                                                        
      placeholder = "Description of the recording device..."                             
    ),                                                                                     
    hr(),                                                                                 
    h5("Export Options"),

    radioButtons('nwb_trace_selection', 'Traces to export:',
      choices = c('All traces' = 'all', 'Selected traces only' = 'selected'),
      selected = 'selected'
    ),
    verbatimTextOutput('tracesInfoNWB'),
    hr(),
    h5("Notes"),
    tags$div(class = "form-group",
      tags$label("Additional Notes"),
      tags$textarea(
        id = "nwb_notes",
        class = "form-control",
        rows = "4",
        placeholder = "Enter any additional experimental notes, observations, or protocol details..."
      )
    ),
    hr(),
    textOutput('nwb_filename_display')
  )
}

create_sidebar_panel <- function() {
  sidebarPanel(
    width = 3,
    fileInput('abfFiles', 'Upload ABF Files', multiple = TRUE, accept = '.abf'),
    tabsetPanel(
      id = "settingsTabs",
      create_main_settings_ui(),
      create_trace_selection_ui(),
      create_graph_settings_ui(),
      create_nwb_settings_ui()
    ),
    hr(),
    actionButton('loadData', 'Load ABF Data', class = 'btn-primary'),
    actionButton('runAnalysis', 'Run Analysis', class = 'btn-primary'),
    hr(),
    downloadButton('downloadExcel', 'Download Excel'),
    downloadButton('downloadRData', 'Download RData'),
    downloadButton('downloadNWB', 'Download NWB'),
    actionButton('clearAll', 'Clear All', class = 'btn-default')
  )
}

create_main_panel <- function() {
  mainPanel(
    width = 9,
    tabsetPanel(
      id = "mainTabs",
      tabPanel("Summary",
        h4("Data Summary"),
        verbatimTextOutput('summaryText'),
        hr(),
        tableOutput('summaryTable')
      ),
      tabPanel("Time Series Plot",
        plotOutput('timeSeriesPlot', height = '800px')
      ),
      tabPanel("Single Examples",
        fluidRow(column(12, downloadButton('downloadSVG', 'Download SVG Plots', class = 'btn-default'))),
        hr(),
        uiOutput('examplePlotsUI')
      ),
      tabPanel("Review Traces",
        fluidRow(
          column(12,
            h4("Review Individual Traces"),
            helpText("Review and accept/reject individual traces for each level")
          )
        ),
        hr(),
        fluidRow(
          column(4,
            selectInput('reviewLevelSelect', 'Select Level to Review:', choices = NULL),
            actionButton('startReview', 'Start Review', class = 'btn-primary'),
            hr(),
            verbatimTextOutput('reviewProgress')
          ),
          column(8,
            plotOutput('reviewPlot', height = '500px'),
            hr(),
            fluidRow(
              column(6, actionButton('acceptTrace', 'Accept', class = 'btn-success', style = 'width: 100%;')),
              column(6, actionButton('rejectTrace', 'Reject', class = 'btn-danger', style = 'width: 100%;'))
            )
          )
        )
      ),
      tabPanel("Data Export",
        h4("Available Data"),
        verbatimTextOutput('exportInfo'),
        hr(),
        h5("Raw Data Preview"),
        tableOutput('rawDataPreview')
      )
    )
  )
}


# DATA PROCESSING HELPER FUNCTIONS


load_abf_files <- function(paths, names) {
  if (length(paths) == 1) {
    readABF(paths[1])
  } else {
    old_wd <- getwd()
    temp_dir <- dirname(paths[1])
    setwd(temp_dir)
    for (i in seq_along(paths)) {
      file.copy(paths[i], names[i], overwrite = TRUE)
    }
    abf_data <- readABFs(names)
    setwd(old_wd)
    abf_data
  }
}

validate_channels <- function(psc_ch, hp_ch, n_channels) {
  if (is.na(psc_ch) || psc_ch < 1 || psc_ch > n_channels) {
    stop(paste("PSC channel", psc_ch, "is invalid. Available channels: 1 to", n_channels))
  }
  if (is.na(hp_ch) || hp_ch < 1 || hp_ch > n_channels) {
    stop(paste("HP channel", hp_ch, "is invalid. Available channels: 1 to", n_channels))
  }
  TRUE
}

extract_channel_data <- function(abf_dataset, channel) {
  data_matrix <- sapply(1:length(abf_dataset$data), function(ii) {
    abf_dataset$data[[ii]][, channel]
  })
  colnames(data_matrix) <- seq(ncol(data_matrix))
  rownames(data_matrix) <- seq(nrow(data_matrix))
  data_matrix
}

detect_stimulation_from_ttl <- function(abf_dataset, ttl_ch, dt, threshold = 0.5) {
  tryCatch({
    TTL_pulse <- sapply(1:length(abf_dataset$data), function(ii) {
      abf_dataset$data[[ii]][, ttl_ch]
    })
    pulse_on <- which(TTL_pulse[, 1] > threshold)
    if (length(pulse_on) > 0) {
      idx2 <- pulse_on[1]
      stim_time <- idx2 * dt - dt
      return(list(success = TRUE, time = stim_time))
    } else {
      return(list(success = FALSE, message = "No TTL pulse detected"))
    }
  }, error = function(e) {
    return(list(success = FALSE, message = e$message))
  })
}

process_traces <- function(data_matrix, dt, stimulation_time, baseline, smooth = 5) {
  out <- lapply(1:ncol(data_matrix), function(ii) {
    peak.fun2(data_matrix[, ii], dt = dt, stimulation_time = stimulation_time, 
             baseline = baseline, smooth = smooth)
  })
  list(
    peak = sapply(out, function(x) x$peak),
    charge = sapply(out, function(x) x$charge),
    response = sapply(out, function(x) x$response),
    baseline = sapply(out, function(x) x$baseline)
  )
}

create_summary_table <- function(Apeak, charge, holding_potential, holding_current) {
  data.frame(
    'time_minutes' = 1:length(Apeak),
    'holding_potential_mV' = holding_potential,
    'holding_current_pA' = holding_current,
    'peak_amplitude_pA' = -Apeak,
    'charge_transfer_pC' = -charge,
    check.names = FALSE
  )
}


# TRACE PARSING AND AVERAGING FUNCTIONS


parse_trace_string <- function(trace_str) {
  if (is.null(trace_str) || nchar(trace_str) == 0) return(integer(0))
  tryCatch({
    if (grepl(':', trace_str)) {
      parts <- strsplit(trace_str, ':')[[1]]
      start_val <- as.integer(parts[1])
      end_val <- as.integer(parts[2])
      if (!is.na(start_val) && !is.na(end_val)) return(start_val:end_val)
    } else {
      vals <- as.integer(strsplit(trace_str, ',')[[1]])
      return(vals[!is.na(vals)])
    }
    integer(0)
  }, error = function(e) integer(0))
}

parse_all_traces <- function(levels, input) {
  traces_list <- list()
  for (i in seq_along(levels)) {
    input_id <- paste0('traces_level_', i)
    traces_list[[i]] <- parse_trace_string(input[[input_id]])
  }
  traces_list
}

calculate_average_traces <- function(I_data2, traces2average, split_include, dt) {
  result <- lapply(seq_along(split_include), function(iii) {
    mask <- split_include[[iii]]
    cols <- seq_along(mask)
    accepted <- cols[mask == 1]
    if (length(accepted) > 0) {
      I_data2[, traces2average[[iii]], drop = FALSE][, accepted, drop = FALSE]
    } else {
      NULL
    }
  })
  non_empty <- !sapply(result, is.null)
  if (!any(non_empty)) return(NULL)
  avg_mat <- sapply(result[non_empty], function(mat) rowMeans(mat))
  if (is.null(dim(avg_mat))) avg_mat <- matrix(avg_mat, ncol = 1)
  avg_mat
}

auto_calculate_limits <- function(single_examples, levels) {
  time <- single_examples[, 'time']
  xlim_common <- c(min(time), max(time))
  all_y <- c()
  for (level in levels) {
    if (level %in% colnames(single_examples)) {
      all_y <- c(all_y, single_examples[, level])
    }
  }
  max_abs_y <- max(abs(all_y), na.rm = TRUE)
  ylim_common <- c(-max_abs_y * 1.1, 5)
  ylim_common[1] <- floor(ylim_common[1] / 10) * 10
  xlim_common[2] <- ceiling(xlim_common[2] / 100) * 100
  list(xlim = xlim_common, ylim = ylim_common)
}


# PLOTTING HELPER FUNCTIONS


plot_peak_amplitude_panel <- function(summary, input, traces2average, split_include) {
  tmax <- max(summary$time_minutes)
  xlim <- c(0, 5 * ceiling(tmax / 5))
  ylim1 <- c(0, ceiling(max(summary$peak_amplitude_pA, na.rm = TRUE) / 100) * 100)
  if (ylim1[2] == 0) ylim1[2] <- 100
  
  par(mar = c(0.5, 6, 3, 1), mgp = c(3.5, 0.7, 0))
  plot(summary$time_minutes, summary$peak_amplitude_pA, type = 'n', xlim = xlim, ylim = ylim1,
       xlab = '', ylab = '', main = '', axes = FALSE, xaxs = 'r', yaxs = 'r', cex.main = input$cex_labels)
  axis(2, at = seq(0, ylim1[2], by = input$ymajor_tick_amp), las = 1, tcl = -0.3, 
       cex.axis = input$cex_axis, lwd = input$lwd_graph)
  mtext(input$ylab, side = 2, line = 4, cex = input$cex_labels)
  points(summary$time_minutes, summary$peak_amplitude_pA, pch = 16, col = input$color, cex = input$cex_points)
  
  if (length(traces2average) > 0) {
    for (iii in seq_along(split_include)) {
      mask <- split_include[[iii]]
      cols <- which(mask == 1)
      if (length(cols) > 0 && length(traces2average[[iii]]) > 0) {
        exY <- summary$peak_amplitude_pA[traces2average[[iii]]][cols]
        exX <- summary$time_minutes[traces2average[[iii]]][cols]
        points(exX, exY, col = 'darkgrey', pch = 16, cex = input$cex_points * 1.2)
      }
    }
  }
  
  if (!is.null(input$drugApplication) && nchar(input$drugApplication) > 0) {
    drug_times <- as.numeric(strsplit(input$drugApplication, ',')[[1]])
    for (i in seq_along(drug_times)) {
      usr <- par('usr')
      rect(xleft = drug_times[i], ybottom = usr[4] - (usr[4] - usr[3]) * 0.05, xright = tmax, ytop = usr[4],
           col = if (i %% 2 == 1) rgb(0.5, 0.5, 0.5, 0.3) else rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    }
  }
}

plot_holding_current_panel <- function(summary, input, traces2average, split_include) {
  tmax <- max(summary$time_minutes)
  xlim <- c(0, 5 * ceiling(tmax / 5))
  ylim2 <- c(floor(min(summary$holding_current_pA, na.rm = TRUE) / 100) * 100, 0)
  
  par(mar = c(0.5, 6, 0.5, 1), mgp = c(3.5, 0.7, 0))
  plot(summary$time_minutes, summary$holding_current_pA, type = 'n', xlim = xlim, ylim = ylim2,
       xlab = '', ylab = '', main = '', axes = FALSE, xaxs = 'r', yaxs = 'r')
  axis(2, at = seq(ylim2[1], 0, by = input$ymajor_tick_hc), las = 1, tcl = -0.3, 
       cex.axis = input$cex_axis, lwd = input$lwd_graph)
  mtext('Holding Current (pA)', side = 2, line = 4, cex = input$cex_labels)
  points(summary$time_minutes, summary$holding_current_pA, pch = 16, col = input$color, cex = input$cex_points)
  
  if (length(traces2average) > 0) {
    for (iii in seq_along(split_include)) {
      mask <- split_include[[iii]]
      cols <- which(mask == 1)
      if (length(cols) > 0 && length(traces2average[[iii]]) > 0) {
        exY <- summary$holding_current_pA[traces2average[[iii]]][cols]
        exX <- summary$time_minutes[traces2average[[iii]]][cols]
        points(exX, exY, col = 'darkgrey', pch = 16, cex = input$cex_points * 1.2)
      }
    }
  }
}

plot_holding_potential_panel <- function(summary, input, traces2average, split_include) {
  tmax <- max(summary$time_minutes)
  xlim <- c(0, 5 * ceiling(tmax / 5))
  ylim3 <- range(summary$holding_potential_mV, na.rm = TRUE)
  ylim3 <- c(floor(ylim3[1] / 10) * 10 - 10, ceiling(ylim3[2] / 10) * 10 + 10)
  
  par(mar = c(5, 6, 0.5, 1), mgp = c(3.5, 0.7, 0))
  plot(summary$time_minutes, summary$holding_potential_mV, type = 'n', xlim = xlim, ylim = ylim3,
       xlab = '', ylab = '', main = '', axes = FALSE, xaxs = 'r', yaxs = 'r')
  axis(2, at = seq(ylim3[1], ylim3[2], by = input$ymajor_tick_hp), las = 1, tcl = -0.3, 
       cex.axis = input$cex_axis, lwd = input$lwd_graph)
  mtext('Holding Potential (mV)', side = 2, line = 4, cex = input$cex_labels)
  axis(1, at = seq(xlim[1], xlim[2], by = input$xmajor_tick_hp), tcl = -0.3, 
       cex.axis = input$cex_axis, lwd = input$lwd_graph)
  mtext(input$xlab, side = 1, line = 2.5, cex = input$cex_labels)
  points(summary$time_minutes, summary$holding_potential_mV, pch = 16, col = input$color, cex = input$cex_points)
  
  if (length(traces2average) > 0) {
    for (iii in seq_along(split_include)) {
      mask <- split_include[[iii]]
      cols <- which(mask == 1)
      if (length(cols) > 0 && length(traces2average[[iii]]) > 0) {
        exY <- summary$holding_potential_mV[traces2average[[iii]]][cols]
        exX <- summary$time_minutes[traces2average[[iii]]][cols]
        points(exX, exY, col = 'darkgrey', pch = 16, cex = input$cex_points * 1.2)
      }
    }
  }
}

draw_scale_bars <- function(xbar, ybar, bar_lwd) {
  usr <- par('usr')
  x_range <- usr[1:2]
  y_range <- usr[3:4]
  ybar_start <- y_range[1] + (y_range[2] - y_range[1]) / 20
  x_start <- x_range[2] - xbar - (x_range[2] - x_range[1]) * 0.05
  y_start <- ybar_start
  x_end <- x_start + xbar
  y_end <- y_start + ybar
  segments(x_start, y_start, x_end, y_start, lwd = bar_lwd, col = 'black')
  segments(x_start, y_start, x_start, y_end, lwd = bar_lwd, col = 'black')
}

plot_single_example <- function(x, y, title, xlim, ylim, xbar, ybar, bar_lwd) {
  par(mar = c(2, 2, 3, 2))
  plot(x, y, type = 'l', col = 'darkgrey', lwd = 1.5, xlim = xlim, ylim = ylim,
       xlab = '', ylab = '', main = title, bty = 'n', axes = FALSE, cex.main = 1.3)
  abline(h = 0, lty = 2, col = 'grey')
  draw_scale_bars(xbar, ybar, bar_lwd)
}


# DOWNLOAD HELPER FUNCTIONS


generate_base_filename <- function(abfFiles) {
  if (!is.null(abfFiles)) {
    file_names <- tools::file_path_sans_ext(abfFiles$name)
    if (length(file_names) == 1) {
      return(file_names[1])
    } else {
      return(paste0(file_names[1], '_', file_names[length(file_names)]))
    }
  } else {
    return(paste0('ABF_', Sys.Date()))
  }
}

create_excel_workbook <- function(state, input) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "summary")
  writeData(wb, "summary", state$summary)
  
   addWorksheet(wb, "raw PSC data")
  writeData(wb, "raw PSC data", as.data.frame(state$I_data))
  
  addWorksheet(wb, "baseline corrected")
  writeData(wb, "baseline corrected", as.data.frame(state$I_data2))
  
  if (!is.null(state$single_examples)) {
    # CHANGE: "Single Examples" → "single examples"
    addWorksheet(wb, "single examples")
    writeData(wb, "single examples", as.data.frame(state$single_examples))
  }
  
  if (length(state$traces2average) > 0) {
    levels <- trimws(strsplit(input$levels, ',')[[1]])
    trace_selection <- data.frame(Level = character(), Accepted = character(), 
                                  Rejected = character(), stringsAsFactors = FALSE)
    for (i in seq_along(levels)) {
      if (length(state$traces2average[[i]]) > 0) {
        accepted_str <- rejected_str <- ""
        if (length(state$split_include) >= i) {
          accepted <- state$traces2average[[i]][state$split_include[[i]] == 1]
          rejected <- state$traces2average[[i]][state$split_include[[i]] == 0]
          if (length(accepted) > 0) accepted_str <- paste(accepted, collapse = ", ")
          if (length(rejected) > 0) rejected_str <- paste(rejected, collapse = ", ")
        } else {
          accepted_str <- paste(state$traces2average[[i]], collapse = ", ")
        }
        trace_selection <- rbind(trace_selection, 
          data.frame(Level = levels[i], Accepted = accepted_str, Rejected = rejected_str))
      }
    }
 
    addWorksheet(wb, "trace selection")
    writeData(wb, "trace selection", trace_selection)
  }
  
  metadata <- data.frame(
    Parameter = c('dt (ms)', 'Stimulation Time (ms)', 'Baseline (ms)', 'N Traces'),
    Value = c(state$dt, state$stimulation_time, input$baseline, ncol(state$I_data))
  )

  addWorksheet(wb, "metadata")
  writeData(wb, "metadata", metadata)
  wb
}


# NWB EXPORT HELPER FUNCTIONS


generate_dandi_nwb_filename <- function(subject_id, session_id = NULL, modality = 'icephys') {
  if (missing(subject_id) || !nzchar(subject_id)) {
    stop('Subject ID required')
  }
  filename <- paste0('sub-', subject_id)
  if (!is.null(session_id) && nzchar(session_id)) {
    filename <- paste0(filename, '_ses-', session_id)
  }
  filename <- paste0(filename, '_', modality, '.nwb')
  return(filename)
}


prepare_nwb_metadata <- function(abf_dataset, channel_index = 1, 
                                 electrode_resistance = '3-5 MOhm', baseline = NULL, stimulation_time = NULL,
                                 traces2save = NULL, level_column_mapping = NULL, notes = NULL) {
  metadata <- extract_metadata(abf_dataset)
  unit <- metadata$channelUnits[channel_index]
  gain_headstage <- metadata$header$fInstrumentScaleFactor[channel_index]
  gain_additional <- metadata$header$fTelegraphAdditGain[channel_index]
  gain_total <- gain_headstage * gain_additional * 1000
  filtering <- metadata$sections$ADCsec[[channel_index]]$fSignalLowpassFilter
  resistance_comp <- metadata$sections$ADCsec[[channel_index]]$fTelegraphAccessResistance
  capacitance_comp <- metadata$sections$ADCsec[[channel_index]]$fTelegraphMembraneCap
  rate <- 1 / metadata$samplingIntervalInSec
  
  trace_data <- abf_dataset$data
  all_traces <- lapply(trace_data, function(mat) mat[, channel_index])
  all_traces <- unname(all_traces)
  
  # Convert level_column_mapping to a Python-compatible dict using reticulate
  if (!is.null(level_column_mapping)) {
    # Wrap each value in list() to prevent reticulate from converting single elements to scalars
    level_column_mapping <- lapply(level_column_mapping, function(x) as.list(as.integer(x)))
    level_column_mapping <- reticulate::dict(level_column_mapping)
  }
  
  # Ensure traces2save is an iterable (list) so reticulate doesn't convert length-1 integer vectors to Python int
  if (!is.null(traces2save)) {
    traces2save_py <- as.list(as.integer(traces2save))
  } else {
    traces2save_py <- NULL
  }

  metadata <- list(
    data = all_traces,
    unit = unit,
    gain_headstage = gain_headstage,
    gain_additional = gain_additional,
    gain_total_mV_per_unit = gain_total,
    filtering = filtering,
    electrode_resistance = electrode_resistance,
    resistance_comp = resistance_comp,
    capacitance_comp = capacitance_comp,
    rate = rate,
    samplingIntervalInSec = metadata$samplingIntervalInSec, 
    baseline_duration_ms = baseline,
    stimulation_time_ms = stimulation_time,
    traces2save = traces2save_py,
    level_column_mapping = level_column_mapping,
    session_notes = notes
  )
}


# create_nwb_file <- function(abf_dataset, nwb_filepath, metadata, 
#                             channel_index = 1, traces2save = NULL, 
#                             baseline = NULL, stimulation_time = NULL,
#                             level_column_mapping = NULL, notes = NULL) {
#   nwb_data <- prepare_nwb_metadata(abf_dataset, channel_index, metadata$electrode_resistance, 
#                                     baseline, stimulation_time, traces2save, level_column_mapping, notes)
  
#   reticulate::py_run_string('
# import numpy as np
# from pynwb import NWBFile, NWBHDF5IO, TimeSeries, ProcessingModule
# from pynwb.file import Subject
# from pynwb.icephys import IntracellularElectrode, VoltageClampSeries
# from datetime import datetime, timezone

# def create_nwb(data_dict, path, meta, baseline_duration_ms=None, stimulation_time_ms=None, traces2average=None, notes=None):
#     rate = float(data_dict["rate"])
#     unit = data_dict["unit"]
#     gain_total = float(data_dict["gain_total_mV_per_unit"])
#     gain_headstage = float(data_dict["gain_headstage"])
#     gain_additional = float(data_dict["gain_additional"])
#     filtering = str(data_dict["filtering"]) + " Hz low-pass"
#     resistance = data_dict["electrode_resistance"]
#     all_data = data_dict["data"]
    
#     n_traces = len(all_data)
#     traces = range(n_traces) if data_dict["traces2save"] is None else data_dict["traces2save"]
    
#     # Build session description with notes if provided
#     session_desc = meta.get("session_description", "Whole-cell voltage clamp")
#     if notes:
#         session_desc = session_desc + ". Notes: " + notes

#     # Create NWB file
#     nwbfile = NWBFile(
#         session_description=session_desc,
#         identifier="NWB_" + meta["subject_id"],
#         session_start_time=datetime.now(timezone.utc),
#         experimenter=[meta["experimenter"]],
#         institution=meta["institution"],
#         lab=meta["lab"],
#         experiment_description=meta.get("experiment_description", ""),  
#         keywords=meta.get("keywords", [])                                
#     )
    
#     # Build subject description
#     subject_desc_parts = []
#     if meta.get("cross"):
#         subject_desc_parts.append("Cross: " + meta["cross"])
#     if meta.get("virus"):
#         subject_desc_parts.append("Virus: " + meta["virus"])
#     if meta.get("virus_injection_site"):
#         subject_desc_parts.append("Injection site: " + meta["virus_injection_site"])
#     subject_description = " | ".join(subject_desc_parts) if subject_desc_parts else "Whole-cell patch clamp recording"
    
#     # Add subject
#     subject = Subject(
#         subject_id=meta["subject_id"],
#         species=meta["species"],
#         age=meta["age"],
#         sex=meta["sex"],
#         genotype=meta["genotype"],
#         description=subject_description
#     )
#     nwbfile.subject = subject
    
#     # Add device
#     device_description = meta.get("device_description", "Patch clamp amplifier for electrophysiology recordings")
#     device = nwbfile.create_device(
#         name=meta["device_name"],
#         description=device_description                  
#     )
    
#     # Add electrode with detailed description
#     electrode = IntracellularElectrode(
#         name="elec0",
#         device=device,
#         cell_id=meta["subject_id"],
#         description="Whole-cell patch-clamp electrode | Original unit: {} | Data stored in amperes | Headstage gain: {:.5g} V/{} | Additional gain: {:.2f}× | Total gain: {:.2f} mV/{}".format(
#             unit, gain_headstage, unit, gain_additional, gain_total, unit),
#         filtering=filtering,
#         location=meta["location"],
#         resistance=resistance
#     )
#     nwbfile.add_icephys_electrode(electrode)
    
#     # Add voltage clamp series for each trace
#     for j, i in enumerate(traces):
#         trace = all_data[i]

#         # Build description
#         series_description = f"Voltage clamp recording - sweep {i}"
#         if traces2average is not None:
#             # add level information if available
#             for level_name, trace_indices in traces2average.items():
#                 #Handle single-element case
#                 if not isinstance(trace_indices, (list, tuple)):
#                     trace_indices = [trace_indices]
#                 if j in trace_indices:
#                     series_description = f"Voltage clamp recording - sweep {i} ({level_name} condition)"


#         # Convert data to amperes if needed (NWB 2.1.0 requirement)
#         if unit == "pA":
#             trace_converted = np.array(trace) * 1e-12  # pA to amperes
#             unit_nwb = "amperes"
#         elif unit == "nA":
#             trace_converted = np.array(trace) * 1e-9  # nA to amperes
#             unit_nwb = "amperes"
#         else:
#             trace_converted = trace
#             unit_nwb = unit

#         series = VoltageClampSeries(
#             name="{}".format(j),
#             data=trace_converted,
#             rate=rate,
#             starting_time=0.0,
#             electrode=electrode,
#             gain=gain_total,
#             capacitance_slow=data_dict["capacitance_comp"],
#             resistance_comp_correction=data_dict["resistance_comp"],
#             description=series_description,               
#             stimulus_description=f"Stimulation time: {stimulation_time_ms} ms; Baseline period: {baseline_duration_ms} ms" if stimulation_time_ms and baseline_duration_ms else "",
#             sweep_number=np.uint32(i),
#             unit=unit_nwb
#         )
#         nwbfile.add_acquisition(series)
    
#     # Add processing module for metadata
#     sweep_module = ProcessingModule(
#         name="icephys", 
#         description="Sweep metadata and analysis parameters from ABF"
#     )
    
#     # Add level column mapping if provided
#     if data_dict.get("level_column_mapping") is not None:
#         import json
#         level_mapping_json = json.dumps(data_dict["level_column_mapping"])
#         mapping_series = TimeSeries(
#             name="level_column_mapping",
#             data=np.array([ord(c) for c in level_mapping_json], dtype=np.int64),
#             unit="none",
#             rate=1.0,
#             starting_time=0.0,
#             description="JSON mapping of level names to column indices: {level_name: [col_indices]}"
#         )
#         sweep_module.add_data_interface(mapping_series)

#     # Add traces2save if provided
#     if data_dict.get("traces2save") is not None:
#         traces_series = TimeSeries(
#             name="traces2save",
#             data=np.array(data_dict["traces2save"], dtype=np.int64),
#             unit="none",
#             rate=1.0,
#             starting_time=0.0,
#             description="Original ABF trace indices exported to NWB (0-indexed)"
#         )
#         sweep_module.add_data_interface(traces_series)
    
#     # Add baseline duration if provided
#     if baseline_duration_ms is not None:
#         baseline_ts = TimeSeries(
#             name="baseline_duration",
#             data=[baseline_duration_ms],
#             unit="milliseconds",
#             timestamps=[0.0],
#             description="Duration of baseline period used for analysis"
#         )
#         sweep_module.add_data_interface(baseline_ts)
    
#     # Add stimulation time if provided
#     if stimulation_time_ms is not None:
#         stim_ts = TimeSeries(
#             name="stimulation_time",
#             data=[stimulation_time_ms],
#             unit="milliseconds",
#             timestamps=[0.0],
#             description="Time of stimulus onset relative to sweep start"
#         )
#         sweep_module.add_data_interface(stim_ts)

#     # Add samplingIntervalInSec (original from ABF)
#     if data_dict.get("samplingIntervalInSec") is not None:
#         sampling_ts = TimeSeries(
#             name="samplingIntervalInSec",
#             data=[data_dict["samplingIntervalInSec"]],
#             unit="seconds",
#             timestamps=[0.0],
#             description="Original sampling interval from ABF file"
#         )
#         sweep_module.add_data_interface(sampling_ts)
    
#     # Add traces for averaging if provided
#     if traces2average is not None and len(traces2average) > 0:
#         avg_info = []
#         level_info = []
#         for level_name, trace_list in traces2average.items():
#             # Handle single-element case where R passes scalar instead of list
#             if not isinstance(trace_list, (list, tuple)):
#                 trace_list = [trace_list]
#             if trace_list:
#                 level_idx = list(traces2average.keys()).index(level_name)
#                 for trace_idx in trace_list:
#                     avg_info.append(trace_idx)
#                     level_info.append(level_idx)
        
#         if len(avg_info) > 0:
#             traces_avg_ts = TimeSeries(
#                 name="traces_for_averaging",
#                 data=np.array(avg_info, dtype=np.int64),
#                 unit="trace_index",
#                 timestamps=np.array(level_info, dtype=np.float64),
#                 description="Trace indices for averaging. Timestamps indicate level number."
#             )
#             sweep_module.add_data_interface(traces_avg_ts)
    
#     nwbfile.add_processing_module(sweep_module)
    
#     # Write to file
#     with NWBHDF5IO(path, "w") as io:
#         io.write(nwbfile)
#     return True
# ')
  
#   py_meta <- list(
#     subject_id = metadata$subject_id,
#     species = metadata$species,
#     age = metadata$age,
#     sex = metadata$sex,
#     genotype = metadata$genotype,
#     cross = metadata$cross,
#     virus = metadata$virus,
#     virus_injection_site = metadata$virus_injection_site,
#     experimenter = metadata$experimenter,
#     institution = metadata$institution,
#     lab = metadata$lab,
#     session_description = metadata$session_description,
#     experiment_description = metadata$experiment_description,
#     keywords = metadata$keywords,
#     device_name = metadata$device_name,
#     device_description = metadata$device_description,
#     location = metadata$location
#   )
  
#   reticulate::py$create_nwb(nwb_data, nwb_filepath, py_meta,
#                            baseline_duration_ms = baseline,
#                            stimulation_time_ms = stimulation_time,
#                            traces2average = nwb_data$level_column_mapping,
#                            notes = notes)
# }

create_nwb_file <- function(abf_dataset, nwb_filepath, metadata, 
                            channel_index = 1, traces2save = NULL, 
                            baseline = NULL, stimulation_time = NULL,
                            level_column_mapping = NULL, notes = NULL) {
  nwb_data <- prepare_nwb_metadata(abf_dataset, channel_index, metadata$electrode_resistance, 
                                    baseline, stimulation_time, traces2save, level_column_mapping, notes)
  
  reticulate::py_run_string('
import numpy as np
from pynwb import NWBFile, NWBHDF5IO, TimeSeries, ProcessingModule
from pynwb.file import Subject
from pynwb.icephys import IntracellularElectrode, VoltageClampSeries
from datetime import datetime, timezone

def create_nwb(data_dict, path, meta, baseline_duration_ms=None, stimulation_time_ms=None, traces2average=None, notes=None):
    rate = float(data_dict["rate"])
    unit = data_dict["unit"]
    gain_total = float(data_dict["gain_total_mV_per_unit"])
    gain_headstage = float(data_dict["gain_headstage"])
    gain_additional = float(data_dict["gain_additional"])
    filtering = str(data_dict["filtering"]) + " Hz low-pass"
    resistance = data_dict["electrode_resistance"]
    all_data = data_dict["data"]
    
    n_traces = len(all_data)
    traces = range(n_traces) if data_dict["traces2save"] is None else data_dict["traces2save"]
    
    session_desc = meta.get("session_description", "Whole-cell voltage clamp")
    if notes:
        session_desc = session_desc + ". Notes: " + notes
    if meta.get("recording_duration"):
      session_desc = session_desc + " Original recording duration: " + str(round(meta["recording_duration"], 1)) + " seconds."

    if meta.get("session_start_time"):
        from dateutil import tz
        session_start_str = str(meta["session_start_time"]).strip()
        if len(session_start_str) >= 19:
            session_start_str = session_start_str[:19]
        session_start_naive = datetime.strptime(session_start_str, "%Y-%m-%d %H:%M:%S")
        # Use timezone from metadata, default to UTC
        tz_str = meta.get("timezone", "UTC")
        local_tz = tz.gettz(tz_str)
        session_start = session_start_naive.replace(tzinfo=local_tz)
    else:
        session_start = datetime.now(timezone.utc)

    nwbfile = NWBFile(
        session_description=session_desc,
        identifier="NWB_" + meta["subject_id"],
        session_start_time=session_start,
        experimenter=[meta["experimenter"]],
        institution=meta["institution"],
        lab=meta["lab"],
        experiment_description=meta.get("experiment_description", ""),
        keywords=meta.get("keywords", [])
    )
    
    subject_desc_parts = []
    if meta.get("cross"):
        subject_desc_parts.append("Cross: " + meta["cross"])
    if meta.get("virus"):
        subject_desc_parts.append("Virus: " + meta["virus"])
    if meta.get("virus_injection_site"):
        subject_desc_parts.append("Injection site: " + meta["virus_injection_site"])
    subject_description = " | ".join(subject_desc_parts) if subject_desc_parts else "Whole-cell patch clamp recording"
    
    subject = Subject(
        subject_id=meta["subject_id"],
        species=meta["species"],
        age=meta["age"],
        sex=meta["sex"],
        genotype=meta["genotype"],
        description=subject_description
    )
    nwbfile.subject = subject
    
    device_description = meta.get("device_description", "Patch clamp amplifier for electrophysiology recordings")
    device = nwbfile.create_device(
        name=meta["device_name"],
        description=device_description
    )
    
    electrode = IntracellularElectrode(
        name="elec0",
        device=device,
        cell_id=meta["subject_id"],
        description="Whole-cell patch-clamp electrode | Original unit: {} | Data stored in amperes | Headstage gain: {:.5g} V/{} | Additional gain: {:.2f}x | Total gain: {:.2f} mV/{}".format(
            unit, gain_headstage, unit, gain_additional, gain_total, unit),
        filtering=filtering,
        location=meta["location"],
        resistance=resistance
    )
    nwbfile.add_icephys_electrode(electrode)
    
    for j, i in enumerate(traces):
        trace = all_data[i]

        series_description = f"Voltage clamp recording - sweep {i}"
        if traces2average is not None:
            for level_name, trace_indices in traces2average.items():
                if not isinstance(trace_indices, (list, tuple)):
                    trace_indices = [trace_indices]
                if j in trace_indices:
                    series_description = f"Voltage clamp recording - sweep {i} ({level_name} condition)"

        if unit == "pA":
            trace_converted = np.array(trace) * 1e-12
            unit_nwb = "amperes"
        elif unit == "nA":
            trace_converted = np.array(trace) * 1e-9
            unit_nwb = "amperes"
        else:
            trace_converted = trace
            unit_nwb = unit

        series = VoltageClampSeries(
            name="{}".format(j),
            data=trace_converted,
            rate=rate,
            starting_time=0.0,
            electrode=electrode,
            gain=gain_total,
            capacitance_slow=data_dict["capacitance_comp"],
            resistance_comp_correction=data_dict["resistance_comp"],
            description=series_description,
            stimulus_description=f"Stimulation time: {stimulation_time_ms} ms; Baseline period: {baseline_duration_ms} ms" if stimulation_time_ms and baseline_duration_ms else "",
            sweep_number=np.uint32(i),
            unit=unit_nwb
        )
        nwbfile.add_acquisition(series)
    
    sweep_module = ProcessingModule(
        name="icephys",
        description="Sweep metadata and analysis parameters from ABF"
    )
    
    if data_dict.get("level_column_mapping") is not None:
        import json
        level_mapping_json = json.dumps(data_dict["level_column_mapping"])
        mapping_series = TimeSeries(
            name="level_column_mapping",
            data=np.array([ord(c) for c in level_mapping_json], dtype=np.int64),
            unit="none",
            rate=1.0,
            starting_time=0.0,
            description="JSON mapping of level names to column indices: {level_name: [col_indices]}"
        )
        sweep_module.add_data_interface(mapping_series)

    if data_dict.get("traces2save") is not None:
        traces_series = TimeSeries(
            name="traces2save",
            data=np.array(data_dict["traces2save"], dtype=np.int64),
            unit="none",
            rate=1.0,
            starting_time=0.0,
            description="Original ABF trace indices exported to NWB (0-indexed)"
        )
        sweep_module.add_data_interface(traces_series)
    
    if baseline_duration_ms is not None:
        baseline_ts = TimeSeries(
            name="baseline_duration",
            data=[baseline_duration_ms],
            unit="milliseconds",
            timestamps=[0.0],
            description="Duration of baseline period used for analysis"
        )
        sweep_module.add_data_interface(baseline_ts)
    
    if stimulation_time_ms is not None:
        stim_ts = TimeSeries(
            name="stimulation_time",
            data=[stimulation_time_ms],
            unit="milliseconds",
            timestamps=[0.0],
            description="Time of stimulus onset relative to sweep start"
        )
        sweep_module.add_data_interface(stim_ts)

    if data_dict.get("samplingIntervalInSec") is not None:
        sampling_ts = TimeSeries(
            name="samplingIntervalInSec",
            data=[data_dict["samplingIntervalInSec"]],
            unit="seconds",
            timestamps=[0.0],
            description="Original sampling interval from ABF file"
        )
        sweep_module.add_data_interface(sampling_ts)
    
    if traces2average is not None and len(traces2average) > 0:
        avg_info = []
        level_info = []
        for level_name, trace_list in traces2average.items():
            if not isinstance(trace_list, (list, tuple)):
                trace_list = [trace_list]
            if trace_list:
                level_idx = list(traces2average.keys()).index(level_name)
                for trace_idx in trace_list:
                    avg_info.append(trace_idx)
                    level_info.append(level_idx)
        
        if len(avg_info) > 0:
            traces_avg_ts = TimeSeries(
                name="traces_for_averaging",
                data=np.array(avg_info, dtype=np.int64),
                unit="trace_index",
                timestamps=np.array(level_info, dtype=np.float64),
                description="Trace indices for averaging. Timestamps indicate level number."
            )
            sweep_module.add_data_interface(traces_avg_ts)
    
    nwbfile.add_processing_module(sweep_module)
    
    with NWBHDF5IO(path, "w") as io:
        io.write(nwbfile)
    return True
')
  
  py_meta <- list(
    subject_id = metadata$subject_id,
    species = metadata$species,
    age = metadata$age,
    sex = metadata$sex,
    genotype = metadata$genotype,
    cross = metadata$cross,
    virus = metadata$virus,
    virus_injection_site = metadata$virus_injection_site,
    experimenter = metadata$experimenter,
    institution = metadata$institution,
    lab = metadata$lab,
    session_description = metadata$session_description,
    session_start_time = metadata$session_start_time,
    recording_duration = metadata$recording_duration, 
    timezone = metadata$timezone,    
    experiment_description = metadata$experiment_description,
    keywords = metadata$keywords,
    device_name = metadata$device_name,
    device_description = metadata$device_description,
    location = metadata$location
  )
  
  reticulate::py$create_nwb(nwb_data, nwb_filepath, py_meta,
                           baseline_duration_ms = baseline,
                           stimulation_time_ms = stimulation_time,
                           traces2average = nwb_data$level_column_mapping,
                           notes = notes)
}


# state management functions

clear_state <- function(state) {
  state$abf_dataset <- NULL
  state$I_data <- NULL
  state$I_data2 <- NULL
  state$holding_potential <- NULL
  state$holding_current <- NULL
  state$stimulation_time <- NULL
  state$dt <- NULL
  state$summary <- NULL
  state$Apeak <- NULL
  state$charge <- NULL
  state$traces2average <- list()
  state$split_include <- list()
  state$single_examples <- NULL
  state$review_level <- NULL
  state$review_index <- 1
  state$review_active <- FALSE
}

reset_inputs <- function(session) {
  updateTextInput(session, 'levels', value = 'control,drug')
  updateTextInput(session, 'drugApplication', value = '')
  updateNumericInput(session, 'baseline', value = 100)
  updateNumericInput(session, 'stimulation', value = 150)
  updateCheckboxInput(session, 'autoDetectStim', value = TRUE)
  updateNumericInput(session, 'pscChannel', value = 1)
  updateNumericInput(session, 'hpChannel', value = 2)
  updateNumericInput(session, 'ttlChannel', value = 3)
  for (i in 1:10) updateTextInput(session, paste0('traces_level_', i), value = '')
}


# MAIN FUNCTION - analyseABF (Modular Version with NWB)


analyseABF <- function() {
  options(shiny.maxRequestSize = 100*1024^2)
  
  ui <- fluidPage(
    create_dark_mode_css(),
    titlePanel("ABF Analysis"),
    sidebarLayout(create_sidebar_panel(), create_main_panel())
  )
  
  server <- function(input, output, session) {
    state <- reactiveValues(
      abf_dataset = NULL, I_data = NULL, I_data2 = NULL,
      holding_potential = NULL, holding_current = NULL,
      stimulation_time = NULL, dt = NULL, summary = NULL,
      Apeak = NULL, charge = NULL, traces2average = list(),
      split_include = list(), single_examples = NULL,
      review_level = NULL, review_index = 1, review_active = FALSE
    )
    
    observeEvent(input$abfFiles, {
      clear_state(state)
      if (!is.null(input$abfFiles)) {
        first_file <- input$abfFiles$name[1]
        base_id <- sub('\\.abf$', '', first_file)
        subject_id <- paste0('m', base_id)
        updateTextInput(session, 'nwb_subject_id', value = subject_id)
      }
    })
    
    observeEvent(input$clearAll, {
      clear_state(state)
      reset_inputs(session)
      showNotification("All data cleared", type = "message", duration = 2)
    })
    
    observeEvent(input$loadData, {
      req(input$abfFiles)
      withProgress(message = 'Loading ABF files...', value = 0, {
        paths <- input$abfFiles$datapath
        names(paths) <- input$abfFiles$name
        incProgress(0.2, detail = "Reading files...")
        
        tryCatch({
          state$abf_dataset <- load_abf_files(paths, input$abfFiles$name)
          incProgress(0.4, detail = "Extracting data...")
          
          state$dt <- state$abf_dataset$samplingIntervalInSec * 1000
          n_channels <- ncol(state$abf_dataset$data[[1]])
          psc_ch <- as.integer(input$pscChannel)
          hp_ch <- as.integer(input$hpChannel)
          ttl_ch <- as.integer(input$ttlChannel)
          
          validate_channels(psc_ch, hp_ch, n_channels)
          state$I_data <- extract_channel_data(state$abf_dataset, psc_ch)
          state$holding_potential <- extract_channel_data(state$abf_dataset, hp_ch)
          
          incProgress(0.6, detail = "Detecting stimulation...")
          
          if (input$autoDetectStim && !is.na(ttl_ch) && ttl_ch > 0 && ttl_ch <= n_channels) {
            result <- detect_stimulation_from_ttl(state$abf_dataset, ttl_ch, state$dt)
            if (result$success) {
              state$stimulation_time <- result$time
              updateNumericInput(session, 'stimulation', value = round(state$stimulation_time, 2))
              showNotification(paste("Auto-detected stimulation at", round(state$stimulation_time, 2), "ms"), type = "message")
            } else {
              state$stimulation_time <- input$stimulation
              showNotification(result$message, type = "warning")
            }
          } else {
            state$stimulation_time <- input$stimulation
          }
          
          if (is.null(state$stimulation_time) || is.na(state$stimulation_time)) {
            stop("Stimulation time must be specified")
          }
          
          incProgress(0.8, detail = "Processing traces...")
          
          psc_processed <- process_traces(state$I_data, state$dt, state$stimulation_time, input$baseline, smooth = 5)
          state$Apeak <- psc_processed$peak
          state$charge <- psc_processed$charge
          state$I_data2 <- psc_processed$response
          rownames(state$I_data2) <- seq(nrow(state$I_data2))
          colnames(state$I_data2) <- seq(ncol(state$I_data2))
          state$holding_current <- psc_processed$baseline
          
          hp_processed <- process_traces(state$holding_potential, state$dt, state$stimulation_time, input$baseline, smooth = 5)
          state$holding_potential <- hp_processed$baseline
          
          state$summary <- create_summary_table(state$Apeak, state$charge, state$holding_potential, state$holding_current)
          incProgress(1.0, detail = "Done!")
          showNotification("ABF data loaded successfully!", type = "message", duration = 5)
          
        }, error = function(e) {
          showNotification(paste("Error loading ABF:", e$message), type = "error", duration = 10)
        })
      })
    })
    
    observeEvent(input$runAnalysis, {
      req(state$I_data2, state$summary)
      withProgress(message = 'Running analysis...', value = 0, {
        levels <- trimws(strsplit(input$levels, ',')[[1]])
        traces2average <- parse_all_traces(levels, input)
        
        if (length(traces2average) == 0 || all(sapply(traces2average, length) == 0)) {
          showNotification("Please specify traces to average in the Trace Selection tab", type = "error", duration = 5)
          return()
        }
        
        incProgress(0.3, detail = "Averaging traces...")
        split_include <- lapply(traces2average, function(traces) rep(1, length(traces)))
        avg_mat <- calculate_average_traces(state$I_data2, traces2average, split_include, state$dt)
        
        if (is.null(avg_mat)) {
          showNotification("No valid traces to average", type = "error", duration = 5)
          return()
        }
        
        time <- seq(nrow(state$I_data2)) * state$dt - state$dt
        state$single_examples <- cbind(time, avg_mat)
        colnames(state$single_examples) <- c('time', levels)
        state$traces2average <- traces2average
        state$split_include <- split_include
        
        incProgress(1.0, detail = "Done!")
        showNotification("Analysis complete!", type = "message", duration = 3)
      })
    })
    
    observeEvent(state$single_examples, {
      req(state$single_examples)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      limits <- auto_calculate_limits(state$single_examples, levels)
      updateNumericInput(session, 'xlim_min_all', value = round(limits$xlim[1], 1))
      updateNumericInput(session, 'xlim_max_all', value = limits$xlim[2])
      updateNumericInput(session, 'ylim_min_all', value = limits$ylim[1])
      updateNumericInput(session, 'ylim_max_all', value = limits$ylim[2])
    })
    
    observe({
      req(state$traces2average)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      updateSelectInput(session, 'reviewLevelSelect', choices = levels, selected = levels[1])
    })
    
    observeEvent(input$startReview, {
      req(state$traces2average, state$split_include)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      level_idx <- which(levels == input$reviewLevelSelect)
      if (length(level_idx) == 0 || length(state$traces2average[[level_idx]]) == 0) {
        showNotification("No traces selected for this level", type = "warning")
        return()
      }
      state$review_level <- level_idx
      state$review_index <- 1
      state$review_active <- TRUE
      showNotification(paste("Reviewing", input$reviewLevelSelect), type = "message")
    })
    
    observeEvent(input$acceptTrace, {
      req(state$review_active, state$review_level)
      level_idx <- state$review_level
      trace_idx <- state$review_index
      state$split_include[[level_idx]][trace_idx] <- 1
      if (trace_idx < length(state$traces2average[[level_idx]])) {
        state$review_index <- trace_idx + 1
      } else {
        state$review_active <- FALSE
        recalculate_averages()
        showNotification("Review complete! Averages updated.", type = "message", duration = 5)
      }
    })
    
    observeEvent(input$rejectTrace, {
      req(state$review_active, state$review_level)
      level_idx <- state$review_level
      trace_idx <- state$review_index
      state$split_include[[level_idx]][trace_idx] <- 0
      if (trace_idx < length(state$traces2average[[level_idx]])) {
        state$review_index <- trace_idx + 1
      } else {
        state$review_active <- FALSE
        recalculate_averages()
        showNotification("Review complete! Averages updated.", type = "message", duration = 5)
      }
    })
    
    recalculate_averages <- function() {
      req(state$I_data2, state$traces2average, state$split_include)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      avg_mat <- calculate_average_traces(state$I_data2, state$traces2average, state$split_include, state$dt)
      if (is.null(avg_mat)) {
        state$single_examples <- NULL
        return(invisible(NULL))
      }
      time <- seq(nrow(state$I_data2)) * state$dt - state$dt
      state$single_examples <- cbind(time, avg_mat)
      non_empty <- !sapply(state$split_include, function(mask) all(mask == 0))
      colnames(state$single_examples) <- c('time', levels[non_empty])
    }
    
    output$fileInfo <- renderPrint({
      req(state$abf_dataset)
      cat("Channels available:", ncol(state$abf_dataset$data[[1]]), "\n")
      cat("Traces:", length(state$abf_dataset$data), "\n")
      cat("Channel names:", paste(state$abf_dataset$channelNames, collapse = ", "), "\n")
      cat("Channel units:", paste(state$abf_dataset$channelUnits, collapse = ", "), "\n")
    })
    
    output$traceSelectionUI <- renderUI({
      req(state$I_data)
      tagList(helpText(paste("Total traces available:", ncol(state$I_data))), uiOutput('dynamicTraceInputs'))
    })
    
    output$dynamicTraceInputs <- renderUI({
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      tagList(
        h5("Plot Axis Limits (applies to all Single Examples)"),
        fluidRow(
          column(6, numericInput('xlim_min_all', 'X min:', value = NULL, step = 10)),
          column(6, numericInput('xlim_max_all', 'X max:', value = NULL, step = 10))
        ),
        fluidRow(
          column(6, numericInput('ylim_min_all', 'Y min:', value = NULL, step = 10)),
          column(6, numericInput('ylim_max_all', 'Y max:', value = NULL, step = 10))
        ),
        hr(),
        h5("Scale Bar Settings"),
        fluidRow(
          column(6, numericInput('xbar_length', 'X bar:', value = 200, min = 1, step = 10)),
          column(6, numericInput('ybar_length', 'Y bar:', value = 25, min = 1, step = 10))
        ),
        fluidRow(column(6, numericInput('bar_lwd', 'Bar thickness:', value = 2, min = 0.5, max = 5, step = 0.5))),
        hr(),
        h5("Trace Selection"),
        lapply(seq_along(levels), function(i) {
          textInput(paste0('traces_level_', i), paste0('Traces for ', levels[i], ':'), placeholder = 'e.g., 1,2,3,4,5 or 1:5')
        })
      )
    })
    
    output$tracesInfo <- renderPrint({
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      traces_list <- parse_all_traces(levels, input)
      cat("Configured trace groups:\n\n")
      for (i in seq_along(traces_list)) {
        if (length(traces_list[[i]]) > 0) {
          cat(levels[i], ": ", paste(traces_list[[i]], collapse = ", "), "\n")
        }
      }
    })
    
    output$tracesInfoNWB <- renderPrint({
      if (input$nwb_trace_selection == 'all') {
        req(state$I_data)
        cat("Exporting ALL", ncol(state$I_data), "traces\n")
      } else {
        req(state$traces2average, state$split_include)
        levels <- trimws(strsplit(input$levels, ',')[[1]])
        cat("Exporting SELECTED traces:\n\n")
        total <- 0
        for (i in seq_along(state$traces2average)) {
          if (length(state$traces2average[[i]]) > 0) {
            mask <- state$split_include[[i]]
            accepted <- state$traces2average[[i]][mask == 1]
            if (length(accepted) > 0) {
              cat(levels[i], ":", paste(accepted, collapse = ", "), "\n")
              total <- total + length(accepted)
            }
          }
        }
        cat("\nTotal:", total, "traces\n")
      }
    })
    
    output$summaryText <- renderPrint({
      req(state$abf_dataset)
      cat("ABF Dataset Information\n=======================\n\n")
      cat("Sampling interval:", state$dt, "ms\n")
      cat("Number of traces:", ncol(state$I_data), "\n")
      cat("Trace length:", nrow(state$I_data), "samples\n")
      cat("Duration:", round(nrow(state$I_data) * state$dt / 1000, 2), "seconds\n")
      if (!is.null(state$stimulation_time)) cat("Stimulation time:", state$stimulation_time, "ms\n")
      if (!is.null(state$summary)) {
        cat("\nPeak amplitude range:", round(min(state$summary$peak_amplitude_pA, na.rm = TRUE), 2), 
            "to", round(max(state$summary$peak_amplitude_pA, na.rm = TRUE), 2), "pA\n")
      }
    })
    
    output$summaryTable <- renderTable({
      req(state$summary)
      head(state$summary, 20)
    }, striped = TRUE, hover = TRUE, bordered = TRUE)
    
    output$timeSeriesPlot <- renderPlot({
      req(state$summary)
      layout(matrix(1:3, ncol = 1), heights = c(3, 2, 2))
      plot_peak_amplitude_panel(state$summary, input, state$traces2average, state$split_include)
      plot_holding_current_panel(state$summary, input, state$traces2average, state$split_include)
      plot_holding_potential_panel(state$summary, input, state$traces2average, state$split_include)
      layout(1)
    })
    
    output$examplePlotsUI <- renderUI({
      req(state$single_examples)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      plot_outputs <- lapply(seq_along(levels), function(i) plotOutput(paste0('examplePlot_', i), height = '400px'))
      do.call(fluidRow, lapply(plot_outputs, function(p) column(6, p)))
    })
    
    observe({
      req(state$single_examples)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      x <- state$single_examples[, 'time']
      xlim_min <- if (!is.null(input$xlim_min_all) && !is.na(input$xlim_min_all)) input$xlim_min_all else min(x)
      xlim_max <- if (!is.null(input$xlim_max_all) && !is.na(input$xlim_max_all)) input$xlim_max_all else max(x)
      
      if (is.null(input$ylim_min_all) || is.na(input$ylim_min_all) || is.null(input$ylim_max_all) || is.na(input$ylim_max_all)) {
        limits <- auto_calculate_limits(state$single_examples, levels)
        ylim_min <- limits$ylim[1]
        ylim_max <- limits$ylim[2]
      } else {
        ylim_min <- input$ylim_min_all
        ylim_max <- input$ylim_max_all
      }
      
      xlim <- c(xlim_min, xlim_max)
      ylim <- c(ylim_min, ylim_max)
      xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
      ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
      bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
      
      for (i in seq_along(levels)) {
        local({
          my_i <- i
          my_level <- levels[my_i]
          output[[paste0('examplePlot_', my_i)]] <- renderPlot({
            if (my_level %in% colnames(state$single_examples)) {
              y <- state$single_examples[, my_level]
              plot_single_example(x, y, my_level, xlim, ylim, xbar, ybar, bar_lwd)
            }
          })
        })
      }
    })
    
    output$reviewProgress <- renderPrint({
      if (!state$review_active) {
        cat("Click 'Start Review' to begin\n")
        return()
      }
      req(state$review_level, state$traces2average)
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      level_name <- levels[state$review_level]
      total_traces <- length(state$traces2average[[state$review_level]])
      current_trace <- state$traces2average[[state$review_level]][state$review_index]
      accepted <- sum(state$split_include[[state$review_level]] == 1)
      rejected <- sum(state$split_include[[state$review_level]] == 0)
      cat("Reviewing:", level_name, "\n")
      cat("Progress:", state$review_index, "/", total_traces, "\n")
      cat("Current trace:", current_trace, "\n\n")
      cat("Accepted:", accepted, "\n")
      cat("Rejected:", rejected, "\n")
    })
    
    output$reviewPlot <- renderPlot({
      req(state$review_active, state$review_level, state$I_data2, state$single_examples)
      level_idx <- state$review_level
      trace_idx <- state$review_index
      actual_trace_num <- state$traces2average[[level_idx]][trace_idx]
      levels <- trimws(strsplit(input$levels, ',')[[1]])
      level_name <- levels[level_idx]
      time <- state$single_examples[, 'time']
      avg_trace <- state$single_examples[, level_name]
      individual_trace <- state$I_data2[, actual_trace_num]
      
      xlim_min <- if (!is.null(input$xlim_min_all) && !is.na(input$xlim_min_all)) input$xlim_min_all else min(time)
      xlim_max <- if (!is.null(input$xlim_max_all) && !is.na(input$xlim_max_all)) input$xlim_max_all else max(time)
      
      if (is.null(input$ylim_min_all) || is.na(input$ylim_min_all) || is.null(input$ylim_max_all) || is.na(input$ylim_max_all)) {
        ylim_range <- range(c(avg_trace, individual_trace), na.rm = TRUE)
        ylim_min <- ylim_range[1] * 1.1
        ylim_max <- 5
      } else {
        ylim_min <- input$ylim_min_all
        ylim_max <- input$ylim_max_all
      }
      
      xlim <- c(xlim_min, xlim_max)
      ylim <- c(ylim_min, ylim_max)
      xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
      ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
      bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
      
      par(mar = c(2, 2, 3, 2))
      plot(time, individual_trace, type = 'l', col = '#CD5C5C', lwd = 2, ylim = ylim, xlim = xlim,
           xlab = '', ylab = '', main = paste(level_name, "- Trace", actual_trace_num),
           bty = 'n', axes = FALSE, cex.main = 1.3)
      abline(h = 0, lty = 2, col = 'grey70')
      lines(time, avg_trace, col = 'darkgrey', lwd = 3)
      draw_scale_bars(xbar, ybar, bar_lwd)
      legend('bottomleft', legend = c('Current Average', paste('Trace', actual_trace_num)),
             col = c('darkgrey', '#CD5C5C'), lwd = c(3, 2), bty = 'n', cex = 1.1)
    })
    
    output$exportInfo <- renderPrint({
      if (is.null(state$summary)) {
        cat("No data loaded yet.\n")
        return()
      }
      cat("Available data for export:\n\n")
      cat("- Summary table (", nrow(state$summary), " rows)\n", sep = "")
      cat("- Raw PSC data (", ncol(state$I_data), " traces)\n", sep = "")
      cat("- Baseline corrected data\n")
      if (!is.null(state$single_examples)) {
        cat("- Single example traces (", ncol(state$single_examples) - 1, " levels)\n", sep = "")
      }
      cat("\nUse download buttons to export data.")
    })
    
    output$rawDataPreview <- renderTable({
      req(state$summary)
      head(state$summary, 10)
    }, striped = TRUE, hover = TRUE)
    
    output$downloadExcel <- downloadHandler(
      filename = function() paste0(generate_base_filename(input$abfFiles), '_analysis.xlsx'),
      content = function(file) {
        req(state$summary, state$I_data, state$I_data2)
        wb <- create_excel_workbook(state, input)
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    output$downloadRData <- downloadHandler(
      filename = function() paste0('ABF_analysis_', Sys.Date(), '.RData'),
      content = function(file) {
        results <- list(
          summary = state$summary, raw_data = state$I_data, baseline_corrected = state$I_data2,
          single_examples = state$single_examples, traces2average = state$traces2average,
          split_include = state$split_include,
          metadata = list(dt = state$dt, stimulation_time = state$stimulation_time, baseline = input$baseline)
        )
        save(results, file = file)
      }
    )

    output$nwb_filename_display <- renderText({
      subject <- if (!is.null(input$nwb_subject_id) && nzchar(input$nwb_subject_id)) input$nwb_subject_id else "Mouse001"
      session <- if (!is.null(input$nwb_session_id) && nzchar(input$nwb_session_id)) input$nwb_session_id else "s01"
      paste0("Filename: sub-", subject, "_ses-", session, "_icephys.nwb")
    })

    output$downloadNWB <- downloadHandler(
      filename = function() {
        generate_dandi_nwb_filename(
          subject_id = input$nwb_subject_id,
          session_id = input$nwb_session_id,
          modality = 'icephys'
        )
      },
      content = function(file) {
        req(state$abf_dataset, state$I_data)
        
        # Step 1: Initialize variables
        traces2save <- NULL
        levs <- NULL
        level_column_mapping <- NULL  # ← FIXED: Initialize here
        
        # Step 2: Handle trace selection
        if (input$nwb_trace_selection == 'selected') {
          req(state$traces2average, state$split_include)
          all_selected <- c()
          all_levels <- c()
          
          for (i in seq_along(state$traces2average)) {
            if (length(state$traces2average[[i]]) > 0) {
              mask <- state$split_include[[i]]
              accepted_idx <- which(mask == 1)
              if (length(accepted_idx) > 0) {
                all_selected <- c(all_selected, state$traces2average[[i]][accepted_idx])
                all_levels <- c(all_levels, rep(i, length(accepted_idx)))
              }
            }
          }
          
          if (length(all_selected) > 0) {
            traces2save <- as.integer(all_selected - 1)  # Convert to 0-based for Python
            
            # Calculate column indices for each level
            level_names <- trimws(strsplit(input$levels, ',')[[1]])
            level_column_mapping <- list()  # Now safe to create
            col_start <- 0
            for (i in seq_along(level_names)) {
              num_traces_in_level <- sum(all_levels == i)
              if (num_traces_in_level > 0) {
                level_column_mapping[[level_names[i]]] <- seq(col_start, col_start + num_traces_in_level - 1)
                col_start <- col_start + num_traces_in_level
              }
            }
          }
        }
        
        # Step 3: Build NWB metadata

        # Parse keywords from comma-separated input
        keywords_raw <- input$nwb_keywords
        keywords <- if (!is.null(keywords_raw) && nzchar(keywords_raw)) trimws(strsplit(keywords_raw, ',')[[1]]) else NULL
        session_description = if (!is.null(input$nwb_session_description) && nzchar(input$nwb_session_description)) input$nwb_session_description else "Whole-cell voltage clamp recording"

        # ABF recording timestamp
        session_start_time <- time_stamp(state$abf_dataset)

        nwb_metadata <- list(
          subject_id = input$nwb_subject_id,
          species = input$nwb_species,
          age = input$nwb_age,
          sex = input$nwb_sex,
          genotype = input$nwb_genotype,
          cross = input$nwb_cross,
          virus = input$nwb_virus,
          virus_injection_site = input$nwb_virus_injection_site,
          experimenter = input$nwb_experimenter,
          institution = input$nwb_institution,
          lab = input$nwb_lab,
          session_description = session_description,
          session_start_time = session_start_time,
          timezone = "America/Chicago", 
          experiment_description = input$nwb_experiment_description,     
          keywords = keywords,                                           
          device_name = input$nwb_device_name,
          device_description = input$nwb_device_description,             
          electrode_resistance = input$nwb_electrode_resistance,
          location = input$nwb_location
        )
        
        # Step 4: Collect analysis parameters
        baseline <- isolate(input$baseline)
        stimulation_time <- isolate(state$stimulation_time)
        notes <- isolate(input$nwb_notes)
        
        # Step 5: Create NWB file with error handling
        tryCatch({
          withProgress(message = 'Creating NWB file...', {
            incProgress(0.5)
            create_nwb_file(
              state$abf_dataset, 
              file, 
              nwb_metadata,
              as.integer(input$pscChannel), 
              traces2save, 
              baseline, 
              stimulation_time, 
              level_column_mapping,  # ← Now always defined (NULL or list)
              notes
            )
          })
          showNotification("NWB file created successfully!", type = "message", duration = 5)
        }, error = function(e) {
          showNotification(paste("NWB Error:", e$message), type = "error", duration = 10)
        })
      }
    )
    
    output$downloadSVG <- downloadHandler(
      filename = function() paste0(generate_base_filename(input$abfFiles), '_plots.zip'),
      content = function(file) {
        req(state$single_examples)
        levels <- trimws(strsplit(input$levels, ',')[[1]])
        x <- state$single_examples[, 'time']
        xlim_min <- if (!is.null(input$xlim_min_all) && !is.na(input$xlim_min_all)) input$xlim_min_all else min(x)
        xlim_max <- if (!is.null(input$xlim_max_all) && !is.na(input$xlim_max_all)) input$xlim_max_all else max(x)
        
        if (is.null(input$ylim_min_all) || is.na(input$ylim_min_all) || is.null(input$ylim_max_all) || is.na(input$ylim_max_all)) {
          limits <- auto_calculate_limits(state$single_examples, levels)
          ylim_min <- limits$ylim[1]
          ylim_max <- limits$ylim[2]
        } else {
          ylim_min <- input$ylim_min_all
          ylim_max <- input$ylim_max_all
        }
        
        xlim <- c(xlim_min, xlim_max)
        ylim <- c(ylim_min, ylim_max)
        xbar <- if (!is.null(input$xbar_length) && !is.na(input$xbar_length)) input$xbar_length else 50
        ybar <- if (!is.null(input$ybar_length) && !is.na(input$ybar_length)) input$ybar_length else 50
        bar_lwd <- if (!is.null(input$bar_lwd) && !is.na(input$bar_lwd)) input$bar_lwd else 2
        base_name <- generate_base_filename(input$abfFiles)
        temp_dir <- tempdir()
        svg_files <- c()
        
        for (level in levels) {
          if (level %in% colnames(state$single_examples)) {
            y <- state$single_examples[, level]
            svg_filename <- file.path(temp_dir, paste0(base_name, '_', level, '.svg'))
            svg_files <- c(svg_files, svg_filename)
            svg(svg_filename, width = 7, height = 5)
            plot_single_example(x, y, level, xlim, ylim, xbar, ybar, bar_lwd)
            dev.off()
          }
        }
        zip(file, svg_files, flags = '-j')
      }
    )
  }
  
  shinyApp(ui, server)
}

# abf2nwb <- function(abf_dataset, nwb_dataset,
#                     subject_id = 'Mouse001',
#                     species = 'Mus musculus',
#                     age = 'P30D',
#                     sex = 'M',
#                     genotype = 'WT',
#                     cross = NULL,
#                     virus = NULL,
#                     virus_injection_site = NULL,
#                     experimenter = 'LastName, FirstName',
#                     institution = 'YourInstitution',
#                     lab = 'YourLab',
#                     session_description = 'Whole-cell voltage clamp',
#                     experiment_description = NULL,                 
#                     keywords = NULL,                               
#                     device_name = 'Multiclamp 700B Patch clamp amplifier',
#                     device_description = NULL,                     
#                     electrode_resistance = '3-5 MOhm',
#                     location = 'Striatum',
#                     channel_index = 1,
#                     traces2save = NULL,
#                     levs = NULL,
#                     baseline = NULL,
#                     stimulation_time = NULL,
#                     notes = NULL) {
  
#   # Build metadata object
#   metadata <- list(
#     subject_id = subject_id,
#     species = species,
#     age = age,
#     sex = sex,
#     genotype = genotype,
#     cross = cross,
#     virus = virus,
#     virus_injection_site = virus_injection_site,
#     experimenter = experimenter,
#     institution = institution,
#     lab = lab,
#     session_description = session_description,
#     experiment_description = experiment_description,      
#     keywords = keywords,                                  
#     device_name = device_name,
#     device_description = device_description,              
#     location = location,
#     electrode_resistance = electrode_resistance
#   )
  
#   # Call create_nwb_file with all parameters
#   create_nwb_file(
#     abf_dataset = abf_dataset,
#     nwb_filepath = nwb_dataset,
#     metadata = metadata,
#     channel_index = channel_index,
#     traces2save = traces2save,
#     baseline = baseline,              
#     stimulation_time = stimulation_time,  
#     level_column_mapping = levs,      
#     notes = notes                     
#   )
# }

abf2nwb <- function(abf_dataset, 
                    nwb_dataset,
                    subject_id = 'Mouse001',
                    species = 'Mus musculus',
                    age = 'P30D',
                    sex = 'M',
                    genotype = 'WT',
                    cross = NULL,
                    virus = NULL,
                    virus_injection_site = NULL,
                    experimenter = 'LastName, FirstName',
                    institution = 'YourInstitution',
                    lab = 'YourLab',
                    session_description = 'Whole-cell voltage clamp',
                    session_start_time = NULL,
                    recording_duration = NULL,   
                    timezone = 'America/Chicago',                    
                    experiment_description = NULL,                 
                    keywords = NULL,                               
                    device_name = 'Multiclamp 700B Patch clamp amplifier',
                    device_description = NULL,                     
                    electrode_resistance = '3-5 MOhm',
                    location = 'Striatum',
                    channel_index = 1,
                    traces2save = NULL,
                    levs = NULL,
                    baseline = NULL,
                    stimulation_time = NULL,
                    notes = NULL) {
  
  # Build metadata object
  metadata <- list(
    subject_id = subject_id,
    species = species,
    age = age,
    sex = sex,
    genotype = genotype,
    cross = cross,
    virus = virus,
    virus_injection_site = virus_injection_site,
    experimenter = experimenter,
    institution = institution,
    lab = lab,
    session_description = session_description,
    session_start_time = session_start_time,             
    recording_duration = recording_duration,
    timezone = timezone,
    experiment_description = experiment_description,      
    keywords = keywords,                                  
    device_name = device_name,
    device_description = device_description,              
    location = location,
    electrode_resistance = electrode_resistance
  )
  
  # Call create_nwb_file with all parameters
  create_nwb_file(
    abf_dataset = abf_dataset,
    nwb_filepath = nwb_dataset,
    metadata = metadata,
    channel_index = channel_index,
    traces2save = traces2save,
    baseline = baseline,              
    stimulation_time = stimulation_time,  
    level_column_mapping = levs,      
    notes = notes                     
  )
}

