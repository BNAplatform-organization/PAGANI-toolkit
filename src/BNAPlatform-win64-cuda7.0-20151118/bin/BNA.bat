@echo *******************************************************************************
@echo *           Welcome to use ParBNA by NICS, EE, Tsinghua University            *
@echo *                            Brain Network Analysis                           *
@echo *          This function calculates characteristics of brain networks,        *
@echo * Degree, Clustering coefficients (Cp), Graph diameter (Lp), Modular structure*
@echo *******************************************************************************

@echo -Please enter a directory containing your networks (.csr files) 
@echo -All the .csr files in the directory will be processed.
@set /p out=^>^><nul
@set /p directory=
@echo.

@echo. >BNA_execute.bat
@echo -Do you want to calculate Degree [Y, N]?
@set /p out=^>^><nul
@set /p mychoice= 
@echo.
@if %mychoice%==Y goto calcu_degree
@if %mychoice%==y goto calcu_degree
@if %mychoice%==N goto ask_Cp
@if %mychoice%==n goto ask_Cp
@:calcu_degree
@echo .\Degree.exe %directory%>>BNA_execute.bat
@goto ask_Cp

@:ask_Cp
@echo -Do you want to calculate Cp [Y, N]?
@set /p out=^>^><nul
@set /p mychoice= 
@echo.
@if %mychoice%==Y goto calcu_Cp
@if %mychoice%==y goto calcu_Cp
@if %mychoice%==N goto ask_Lp
@if %mychoice%==n goto ask_Lp
@:calcu_Cp
@echo -How many random networks for Cp comparison? Enter 0 if desiring no comparison.
@set /p out=^>^><nul
@set /p number= 
@echo.
@echo .\Cp.exe %directory% %number% >>BNA_execute.bat
@goto ask_Lp

@:ask_Lp
@echo -Do you want to calculate Lp [Y, N]?
@set /p out=^>^><nul
@set /p mychoice= 
@echo.
@if %mychoice%==Y goto calcu_Lp
@if %mychoice%==y goto calcu_Lp
@if %mychoice%==N goto ask_bc
@if %mychoice%==n goto ask_bc
@:calcu_Lp
@echo -How many random networks for Lp comparison? Enter 0 if desiring no comparison.
@set /p out=^>^><nul
@set /p number= 
@echo.
@echo .\CUBFW_Lp.exe %directory% %number% >>BNA_execute.bat
@goto ask_bc

@:ask_bc
@echo -Do you want to calculate Betweenness Centrality [Y, N]?
@set /p out=^>^><nul
@set /p mychoice= 
@echo.
@if %mychoice%==Y goto calcu_bc
@if %mychoice%==y goto calcu_bc
@if %mychoice%==N goto end
@if %mychoice%==n goto end
@:calcu_bc
@echo .\CUBC.exe %directory%>>BNA_execute.bat

@:end
@echo @pause>>BNA_execute.bat
@.\BNA_execute.bat
@pause