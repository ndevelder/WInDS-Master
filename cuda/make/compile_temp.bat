@ECHO OFF
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\Tools\vsvars32.bat"
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v4.1\bin\nvcc.exe" -I "C:\Program Files\MATLAB\R2011b\extern\include" --ptxas-options=-v -arch sm_20 --cuda "C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\source\BiotSavartMex.cu" --output-file "C:\Users\admin\Documents\MATLAB\Nate\WInDS-2-2012\cuda\compile_temp\BiotSavartMex.cpp"
