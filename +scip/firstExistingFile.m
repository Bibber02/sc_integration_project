function filePath = firstExistingFile(candidates)
%FIRSTEXISTINGFILE Return the first existing file, or the canonical path.

if ischar(candidates) || isstring(candidates)
    candidates = cellstr(candidates);
end

for k = 1:numel(candidates)
    candidate = char(candidates{k});
    if isfile(candidate)
        filePath = candidate;
        return;
    end
end

filePath = char(candidates{1});
end
