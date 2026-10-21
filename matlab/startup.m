pth = mfilename('fullpath');
dir = fileparts(pth);
addpath(dir);
addpath([dir '/src'])
if exist([dir '/../FMM3D/matlab'], 'dir')
    addpath([dir '/../FMM3D/matlab']);
end
