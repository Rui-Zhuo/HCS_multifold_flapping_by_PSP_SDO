function cfg = hcs_config()
%HCS_CONFIG Central paths and reproducibility settings for MATLAB scripts.

cfg.random_seed = 20210117;
repository_root = fileparts(mfilename('fullpath'));
cfg.psp_data_root = env_or_default('HCS_PSP_DATA_ROOT', ...
    fullfile(repository_root, 'data', 'PSP'));
cfg.work_root = env_or_default('HCS_WORK_ROOT', ...
    fullfile(repository_root, 'outputs'));
cfg.legacy_work_root = env_or_default('HCS_LEGACY_WORK_ROOT', ...
    fullfile(repository_root, 'outputs', 'legacy'));
cfg.psi_data_root = env_or_default('HCS_PSI_DATA_ROOT', ...
    fullfile(repository_root, 'data', 'PSI'));
cfg.helioviewer_data_root = env_or_default('HCS_HELIOVIEWER_DATA_ROOT', ...
    fullfile(repository_root, 'data', 'HelioViewer'));

if exist('hcs_config_local', 'file') == 2
    local_cfg = hcs_config_local();
    local_fields = fieldnames(local_cfg);
    for i_field = 1:numel(local_fields)
        cfg.(local_fields{i_field}) = local_cfg.(local_fields{i_field});
    end
end
end

function value = env_or_default(variable_name, default_value)
value = getenv(variable_name);
if isempty(value)
    value = default_value;
end
end
