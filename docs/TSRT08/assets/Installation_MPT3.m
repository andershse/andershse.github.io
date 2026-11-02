
%% Installation of tbxmanager with all modules required for MPT3
%

clc;
disp('----------------------------------------------');
disp('Installation of MPT using the Toolbox manager.');
disp('----------------------------------------------');
disp(' ');

% Check MATLAB compatibility
if verLessThan('matlab','9.12')
    error('The updated Toolbox Manager requires MATLAB R2022a or newer.');
end

fprintf(['Choose the installation directory where to install the Toolbox manager.\n',...
      'A folder "tbxmanager" will be created in the specified location.\n',...
      'If you do not specify the folder, the current directory will be used.\n']);

% get the installation folder
default_dir = pwd;
c = uigetdir(pwd);

if isequal(c,0)
    fprintf(['No directory has been provided.\n',...
        'Installing the toolbox manager in the current directory "%s".\n'],default_dir);
    c = default_dir;
end

% create a new directory in that folder
d = fullfile(c,'tbxmanager');

if exist(d,'dir') == 7
    disp('The directory "tbxmanager" already exists. Reusing it.');
else
    disp('Creating the directory "tbxmanager".');
    out = mkdir(d);

    if ~out
        error('Could not create the installation directory "%s".',d);
    end
end

% enter that directory
cd(d);

% remove MPT2 or YALMIP from MATLAB path
disp(' ');
disp('Removing toolboxes that may conflict with MPT from the Matlab path.');

if ~isempty(which('mpt_init'))
    rmpath(genpath(fileparts(which('mpt_init'))));
end

if ~isempty(which('yalmipdemo'))
    rmpath(genpath(fileparts(which('yalmipdemo'))));
end

% download the updated tbxmanager
disp(' ');
disp('Downloading the updated Toolbox manager from the internet.');

websave('tbxmanager.m', ...
    'https://marekwadinger.github.io/tbxmanager/tbxmanager.m');

% activate updated Toolbox Manager
addpath(d,'-begin');
clear tbxmanager
rehash;

disp('Active Toolbox Manager:');
which tbxmanager -all

if ~strcmp(which('tbxmanager'),fullfile(d,'tbxmanager.m'))
    error('MATLAB is not using the updated Toolbox Manager.');
end

% install all required modules
disp(' ');
disp('Installing MPT3 and all required modules.');

tbxmanager install mpt mptdoc --yes

% enable MPT3 and all required modules
disp(' ');
disp('Enabling MPT3 and its dependencies.');

tbxmanager enable mpt mptdoc yalmip sedumi lcp cddmex glpkmex clpmex espresso hysdel oasesmex qpspline

% restore installed toolbox paths
disp(' ');
disp('Restoring toolbox paths.');

tbxmanager restorepath
rehash;

% verify MPT3 availability
if isempty(which('mpt_init'))
    error(['MPT3 initialization file was not found. ', ...
           'Check the installation and enabled packages.']);
end

% create the initialization file to set the path
disp(' ');
disp('Creating the initialization file "startup.m".');

p = which('startup.m');

if isempty(p)
    userDirs = strsplit(userpath,pathsep);
    if ~isempty(userDirs{1})
        p = fullfile(userDirs{1},'startup.m');
    else
        p = fullfile(d,'startup.m');
    end
end

% Avoid duplicating initialization lines
marker = '% MPT3 Toolbox Manager Setup';

if exist(p,'file') == 2
    existingText = fileread(p);
else
    existingText = '';
end

if ~contains(existingText,marker)

    fid = fopen(p,'a');

    if isequal(fid,-1)
        error('Could not modify the initialization file "%s".',p);
    end

    fprintf(fid,'\n%s\n',marker);
    fprintf(fid,'addpath(''%s'', ''-begin'');\n', ...
        strrep(d,'''',''''''));
    fprintf(fid,'tbxmanager restorepath;\n');

    fclose(fid);
end

disp('Startup file has been configured.');

% get back to the original directory
cd(default_dir);

% add path to tbxmanager
disp(' ');
disp('Adding path to Matlab.');
addpath(d,'-begin');

% save path for future sessions
disp(' ');
disp('Saving path for future sessions.');

status = savepath;

if status
    warning(['Could not save the MATLAB path automatically. ', ...
             'The startup.m configuration may still restore the toolboxes.']);
end

disp(' ');
disp('Installation finished.');
disp('Next time you start Matlab the toolboxes will be automatically initialized.');

% initialize MPT
disp(' ');
disp('Initializing MPT3.');

mpt_init;

% test MPT3 installation
disp(' ');
disp('Testing MPT3 installation.');

P = Polyhedron('lb',[-1;-1],'ub',[1;1]);

figure;
plot(P);
grid on;
title('MPT3 Installation Test');

disp('----------------------------------------------');
disp('MPT3 installation and basic test completed.');
disp('----------------------------------------------');
