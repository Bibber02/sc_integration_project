function p = paths()
%PATHS Central project paths for the integration project.

root = fileparts(fileparts(mfilename('fullpath')));

p.root = root;
p.model = fullfile(root, 'model');
p.controllers = fullfile(root, 'controllers');
p.lqr = fullfile(p.controllers, 'lqr');
p.kalman = fullfile(p.controllers, 'kalman');
p.robustV1 = fullfile(p.controllers, 'robust_v1');

p.sensorCalibration = fullfile(root, 'sensor_calibration');
p.systemIdentification = fullfile(root, 'system_identification');
p.sensorNoiseMeasurement = fullfile(p.systemIdentification, 'sensor_noise_measurement');
p.fullSystem = fullfile(p.systemIdentification, 'full_system');
p.fullSystemMeasurementData = fullfile(p.fullSystem, 'measurement_data');
p.fullSystemBlackBox = fullfile(p.fullSystem, 'black_box');
p.fullSystemGreyBox = fullfile(p.fullSystem, 'grey_box');
p.fullSystemGreyBoxViscous = fullfile(p.fullSystemGreyBox, 'viscous');
p.fullSystemGreyBoxStribeck = fullfile(p.fullSystemGreyBox, 'stribeck');

p.passiveLink = fullfile(p.systemIdentification, 'passive_link');
p.passiveLinkMeasurementData = fullfile(p.passiveLink, 'measurement_data');
p.passiveLinkGreyBox = fullfile(p.passiveLink, 'grey_box');
p.passiveLinkModels = fullfile(p.passiveLinkGreyBox, 'models');
p.passiveLinkViscous = fullfile(p.passiveLinkGreyBox, 'passive_link_id_viscous');
p.passiveLinkStribeck = fullfile(p.passiveLinkGreyBox, 'passive_link_id_3step_stribeck');

p.templateRotatingPendulum = fullfile(root, 'template', 'rotating-pendulum');
p.hardwareRotatingPendulum = fullfile(root, 'hardware', 'rotating_pendulum');

p.linearizedPlantUprightTs001 = fullfile(p.model, 'linearized_plant_eq_pi_pi_ts_0.001.mat');
p.fullSystemStribeckResultWithP0 = scip.firstExistingFile({
    fullfile(p.fullSystemGreyBoxStribeck, 'results_with_p0', 'full_system_id_new_data_locked_passive_result.mat')
});
p.passiveLinkViscousResult = scip.firstExistingFile({
    fullfile(p.passiveLinkViscous, 'passive_link_viscous_result.mat')
    fullfile(p.passiveLinkViscous, 'passive_link_viscous_result_v5_no_helpers.mat')
});
p.passiveLinkStribeckResult = scip.firstExistingFile({
    fullfile(p.passiveLinkStribeck, 'passive_link_3step_stribeck_result.mat')
    fullfile(p.passiveLinkStribeck, 'passive_link_3step_stribeck_result_v5_no_helpers.mat')
});
end
