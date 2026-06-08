function signature = kalman_reference_signature(settings)
%KALMAN_REFERENCE_SIGNATURE Cache key for the identified reference model.

signature = struct();
signature.TsData = settings.TsData;
signature.inputSign = settings.inputSign;
signature.amplitudes = settings.amplitudes;
signature.idxReference = settings.idxReference;
signature.usePrbs = settings.usePrbs;
signature.useChirp = settings.useChirp;
signature.model = 'greybox_id_full_stribeck_model';
end

