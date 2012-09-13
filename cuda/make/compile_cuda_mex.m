function [ ] = compile_cuda_mex(input_filename)
%compile_cuda_mex Compiles given .cu file first to .cpp and then to mex
%function
%System dependent changes might to be made based on the locations of nvcc
%and the mex extern folder and CUDA toolkit/lib



%get extension-less name of input file
name = strrep(input_filename, '.cu', '')

% Define location of visual studio environment setup
envsetup = '"C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\Tools\vsvars32.bat"';
nvccloc = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.1\bin\nvcc.exe"';
matlabextern = '"C:\Program Files\MATLAB\R2011b\extern\include"';
cudainput = ['"C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\source\' input_filename '"'];
coutput = ['"C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\compile_temp\' name '.cpp"'];
cudatoolkitinc = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.1\include"';
cudalib = '"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.1\lib\x64"';
mexinput = ['"C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\compile_temp\' name '.cpp"'];
outdir = 'C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\core';

%create batch line to setup environment
set_env_vs = ['call ' envsetup];

%create batch line to call nvcc
compile_cuda = [nvccloc ' -I ' matlabextern ' --ptxas-options=-v -arch sm_20 --cuda ' cudainput ' --output-file ' coutput];

%create initial batch line
start = '@ECHO OFF';

%write batch lines to batch file
fileID = fopen('C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\make\compile_temp.bat','w', 'n', 'US-ASCII');
fprintf(fileID,'%s\r\n',start);
fprintf(fileID,'%s\r\n',set_env_vs);
fprintf(fileID,'%s\r\n',compile_cuda);
fclose(fileID);

%run the batch file
run_bat = system('C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\make\compile_temp.bat');

%feedback on batch file run
if(run_bat == 0)
    display('Successfully ran the batch file');
    %create mex function from compiled code
    eval(['clear functions']);
    eval(['mex -I' cudatoolkitinc ' -L' cudalib ' -lcudart -lcuda ' mexinput ' -outdir ' outdir ' COMPFLAGS="$COMPFLAGS -openmp"  LINKFLAGS="$LINKFLAGS -openmp"']);
end



end





