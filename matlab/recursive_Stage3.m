% ESTIMATES THE HLW model with both correct and incorrct Stage 2 procedures, as well as
% \sigma(g) and \sigma(z) directly.
% ---------------------------------------------------------------------------------------------------
% Reminder: their SSF is
% ---------------------------------------------------------------------------------------------------
% 	Y(t) = A*X(t) + H*Xi(t) + v(t), Var[v(t)] = R.
%  Xi(t) = F*Xi(t-1) + S*e(t),			Var[e(t)] = Q.
% ---------------------------------------------------------------------------------------------------
% MY STATE SPACE FORM IS:
% 		Observed:	Y(t)			= D(t) + M*alpha(t)			+ e(t);		Var[e(t)] = H.
% 		State:		alpha(t)	= C(t) + Phi*alpha(t-1)	+ S*n(t);	Var[n(t)] = Q.
% ===================================================================================================
clear; clc;
TT = (110:1:240);
% TT = 210:1:213;
% TT = ([162:1:180]);
TT = (110);

set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 2 in plots;
%% ADD LOCAL FUNCTION PATH WITHOUT SUBFOLDERS (IE. '_OLD') (only add what is needed)
addpath('local.Functions','utility.Functions');
cntr_ = {'US', 'EA', 'UK', 'CA'};

% MAKE OUTPUT/GRAPHICS DIRECTORIES IF THEY DO NOT EXIST
matlab_output_dir		= '../matlab.output/';
latex_graphics_dir	= '../../graphics/';
latex_table_dir			= '../../table.input/';
% making dirs
if ~exist(matlab_output_dir,	'dir');	mkdir(matlab_output_dir);		end
if ~exist(latex_graphics_dir,	'dir');	mkdir(latex_graphics_dir);	end
if ~exist(latex_table_dir,		'dir');	mkdir(latex_table_dir);			end

% DEFINDE WHICH VINTAGE OF DATA TO BE USE IN ESTIMATION. THESE ARE STORED IN DIFFERENT DIRECTORIES, 
% DATA_DIR_INPUT	= '../data/R.data.for.estimation.2021.Jan.29/';		% DATA ENDS IN Q4-2019
DATA_DIR_INPUT	= '../data/R.data.for.estimation.2020.May.28/';		% DATA ENDS IN Q4-2019
% DATA_DIR_INPUT	= '../data/R.data.for.estimation.2020.Oct.5/';	% DATA ENDS IN Q2-2020

% for CI = 1:4
CI = 1; 
	% DEFINE COUNTRY, IE 
	COUNTRY		= cntr_{CI};
	% DEFINE SAMPLE-END. NOTE IF DATA IN CSV/MAT FILE WHICH IS READ IS LESS THAN SMPL_END, DATA ENDS THERE.
	SMPL_END	= 'Q4-2019'; 
	% SET VARIOUS PLOTTING AND DATA READING/WRITING 
	CSV_READ = 0;				% SET TO 1 TO READ NEW DATA FROM CSV FILE, OTHERWISE LOAD THE .MAT CONVERTED FILE
	% KEEP A LOG FILE USING THE DIARY FUNCTION
	DIARY_ON = 0;

	% WRITE RESULTS TO LOG FILE
	if DIARY_ON
		% save as textfile
    %   diary_file_name = ['./_diary_dir/' char(cntr_(CI)) '_results.txt'];
		% save to _output folder as well
		diary_file_name = [matlab_output_dir char(cntr_(CI)) '_results.txt'];
		if exist(diary_file_name,'file'); delete(diary_file_name); end
		diary(diary_file_name)
	end
	
	%% OPTIMISATION SETTINGS FOR FMINUNC AND FMINCON
	% ***************************************************************************************************
	% NOTE: STAGE 1 ESTIMATES SENSITIVE TO INITIAL VALUES AND OPTIMISATION SETTINGS
	% ---------------------------------------------------------------------------------------------
	% FMINUNC HAS FOLLOWING ALGORITHMS: 
	% ---------------------------------------------------------------------------------------------
	% (1) 'quasi-newton' (default); (2) 'trust-region'.
	% ---------------------------------------------------------------------------------------------
	% Use optimoptions to set the Algorithm option at the command line.Recommendations: If your 
	% objective function includes a gradient, use 'Algorithm' = 'trust-region', and set 
	% the SpecifyObjectiveGradient 
	% option to true. Otherwise, use 'Algorithm' = 'quasi-newton'.
	% ---------------------------------------------------------------------------------------------
	options_unc = optimoptions(@fminunc, ...
								'Algorithm'								,	'quasi-newton' , ... 
								'Display'									, 'none'	, ...
								'MaxIterations'						, 5e5			, ... 
								'MaxFunctionEvaluations'	, 5e5			, ... 
								'OptimalityTolerance'			,	1e-8		, ...
								'FunctionTolerance'				, 1e-8		, ...
								'StepTolerance'						, 1e-8		, ...
								'FiniteDifferenceType'		, 'central');
	% 						'HessianApproximation'		, 'bfgs'	, ...						
	% 						'Algorithm'		, 'trust-region', ...
	% 						'Algorithm'		, 'quasi-newton', ...
	% ---------------------------------------------------------------------------------------------
	% FMINCON HAS FOLLOWING ALGORITHMS: 
	% ---------------------------------------------------------------------------------------------
	% (1) 'interior-point' (default); (2) 'trust-region-reflective'; (3) 'sqp'; (4) 'sqp-legacy';
	% (5) 'active-set'
	% ---------------------------------------------------------------------------------------------
	% Use the 'interior-point' algorithm first.
	% For help if minimization fails, see When the Solver Fails or When the Solver Might Succeed. 
	% Try 'sqp' next, and 'active-set' last. Use 'trust-region-reflective' when applicable.	
	% ---------------------------------------------------------------------------------------------
	% Set Algorithm to 'active-set' to get their estimates, 
	% changing this to 'interior-point', to get the higher loglike value
	% ---------------------------------------------------------------------------------------------
	options_con = optimoptions(@fmincon, ... 
								'Algorithm'								,	'interior-point'		, ... 
								'Display'									, 'none'	, ...
								'MaxIterations'						, 5e5			, ... 
								'MaxFunctionEvaluations'	, 5e5			, ... 
								'HessianApproximation'		, 'bfgs'	, ...
								'OptimalityTolerance'			,	1e-8		, ...
								'FunctionTolerance'				, 1e-8		, ...
								'StepTolerance'						, 1e-8		, ...
								'FiniteDifferenceType'		, 'central');
	% 						'Algorithm'		, 'trust-region-dogleg', ...
	% 						'Algorithm'		, 'active-set', ...												
	% 						'Algorithm'		, 'quasi-newton', ...
	% ***************************************************************************************************
	% PARAMETER NAMES
	par_names3 = {'a_y1'    ; ...
								'a_y2'    ; ...
								'a_r'     ; ...
								'b_pi'    ; ...
								'b_y'     ; ...
								'sigma_y~'; ...
								'sigma_pi'; ...
								'sigma_y*'; ...
								'sigma_g' ; ...
								'sigma_z' ; ...
								'Log-Like'};
% 							; ...
% 								'Lambda_g'; ...
% 								'Lambda_z'};
	% ***************************************************************************************************

	%% LOAD THE HWL INPUT DATA TABLE
	% ---------------------------------------------------------------------------------------------------
	if strcmp(COUNTRY,'EA') 
		hlw_read_sample = '1972Q1.to.2019Q4';	
	else 
		hlw_read_sample = '1961Q1.to.2019Q4'; 
	end
	
	DATA_DIR_NAME	= [COUNTRY '.data'];
	%	READ IN DATA FROM MATLAB OR CSV FILE
	data = read_data_R_csv(DATA_DIR_NAME, DATA_DIR_INPUT, CSV_READ);
	% head2tail(data)
	% MAKE ANONYMOUS FUNCTION TO GET SIGMA/LAMBDA
	S3g	= @(theta,L_g) abs(theta(8)*L_g);
	S2g	= @(theta,L_g) abs(theta(10)*L_g);
	L2g	= @(theta) theta(11)/theta(10);
	L3g	= @(theta) theta(9)/theta(8);
	L3z	= @(theta) abs( theta(10)*theta(3)/theta(6) );
	S3z	= @(theta,L_z) abs(theta(6)*L_z/theta(3));
	% now for reading of some data
	Sa00 = @(x) (['Stage' num2str(x) '.xi00.']	);		% 'x00.';
	SP00 = @(x) (['Stage' num2str(x) '.P00.']	);			% 'P00.';
	Sstr = @(x) (['Stage' num2str(x) '.theta.']);			% 'theta.';
	S2Lz = @(x) (['Stage' num2str(x) '.Lambda.z.']);	% 'Lambda.z.';
	% HLW R-code data output of parameters, some data etc, which is used as an input.
	HLW_INPUT_DATA_DIR = '../data/R.HLW.results/';

	% CHANGE SAMPLE SIZE IF NEEDED (use timerange function to shorten sample).
	data_full = data(timerange('Q1-1960',SMPL_END,'closed'),:);
	% data = data(timerange('Q1-1960','Q4-2007','closed'),:); % does not matter if beg-sample is less than in data for EA.
	
	%% make storage for stage 3 MLE parameter estimates
	S3MLE = [];
 	for i = 1:length(TT)
	% TRIM THE DATA TO RIGHT SIZE
 	data = data_full(1:TT(i),:);
% 	data = data_full;

% % % 
% % % 	% ---------------------------------------------------------------------------------------------------
% % % 	% HLW: LOAD PARAMETER ESTIMATES FROM THEIR R FILES.
% % % 	% ---------------------------------------------------------------------------------------------------
% % % 	% function handle to load the HLW R output paramters etc
% % % 	getS = @(x) ([HLW_INPUT_DATA_DIR COUNTRY '/' x '.csv']);
% % % 	HLW_Stage3_R_File	= xlsread( getS([Sstr(3) hlw_read_sample]) );
% % % 	% THESE ARE THE INTIAL VALUES USED BY HLW (FIRST COLUMN OF HLW2017_THETA)
% % % 	initVals3_R_File	= HLW_Stage3_R_File(1:8,1);
% % % 	% these are their estimates
% % % 	HLW_theta_R_File	= HLW_Stage3_R_File(1:8,2);
% % % 
% % % 	% SET OTHER OPTIONAL PARAMETERS NEEDED IN PROCEDURE (MW,EW,QLR, Lambda_z Col(1) MUE.stat Col(2))
% % % 	% HLW_S2_Lz_	= dlmread([HLW_input_data_dir 'Stage2.Lambda_z.' SMPL_END '.csv'],',',1,1);
% % % 	Lambda_g = HLW_Stage3_R_File(end-1,2);
% % % 	Lambda_z = HLW_Stage3_R_File(end  ,2);
% % % 
% % % 	% Stage 3 state vector prior (initvals) to be parsed to the function as other paramters
% % % 	a00_S3_R_File = dlmread( getS([Sa00(3) hlw_read_sample]) ,',',0,0);
% % % 	P00_S3_R_File = dlmread( getS([SP00(3) hlw_read_sample]) ,',',0,0);
	% SIMPLE 0.2*I AS IN THEIR INITIALISATIONS
% 	other_pars_S3.P00 = 0.2*eye(size(other_pars_S3.a00,1));
	
% 	other_pars_S3.Lambda_g	= Lambda_g;			%%% some reference values 0.0538690377714457
% 	other_pars_S3.Lambda_z	= Lambda_z;			%%% some reference values 0.0302172212554928

	%% ARANANGE DATA INTO APPROPRIATE Y AND X
	% ---------------------------------------------------------------------------------------------------
	% prepare data for the matrices
	lnGDP	= 100*data.gdp_log;			        % ln(GDP) times 100: NOTE GDP already logged in HLW's R-File
	INFL	=			data.inflation;				    % INFLATION SERIES UNTRANSFORMED
	RR		=			data.real_rate;						% REAL RATE: 

	% defined lagged variables to be used in stage 3 SSF
	GPD_2_lags	= mlag(lnGDP,2);	        % [y_{t-1} y_{t-2}]
	RR_2_lags		= mlag(RR,2);			        % [r_{t-1} r_{t-2}]
	INFL_4_lags = mlag(INFL,4);		        % [\pi_{t-1} \pi_{t-2} \pi_{t-3} \pi_{t-4}]

	% MEASUREMENT MATRIX (2xT) TO BE PASSED TO KF ROUTINE. NOTE THE DATA TRIMMED TO (5:END)
	% (SO FROM 1961:Q1 ONWARDS)
	YY	= [lnGDP INFL]; 
	Y		=	YY(5:end,:);
	% EXOGENEOUS VECTOR X NEEDED FOR DT = A*X (7xT) 
	% [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
	XX	= [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ];
	X		= XX(5:end,:);   
	% sample size
	T		= length(Y);

	% PRINT DATES OF TIME PERIOD CONSIDERED (5:end) is 1961:Q1.
	Dates	= datenum(data.Time(5:end));
	sep('=')
	fprintf( '		Country is %s: ', cntr_{CI}) 
	fprintf( ' Sample period is: %s', datestr(Dates(1),'yyyy:qq'));
	fprintf( ' to %s ', datestr(Dates(end),'yyyy:qq'));
	fprintf( '| Vintage of Data is %s: \n', DATA_DIR_INPUT(31:end-1) );
	sep('=')

	%% SET INITIAL VALUES OF PARAMETERS FOR STAGE 3 MODELS (based on sample going back to 1960:Q1)
	% ---------------------------------------------------------------------------------------------------
	% run HP Filter on Y to get initial crude trend and cycle decomposition
	[HP_cycle, HP_trend] = hp_filter(YY(:,1), 36000); % Lambda = 36000 is in HLW paper
	% HP based trend growth
	delta_HP_trend = delta(HP_trend);
	inflt = YY(:,2) ; 
	% RUN OLS ON HP YCYCL AND RT_1 TO GET INTIAL ESTIMATES A_ AND B_ AND THEIR STANDARD ERRORS
	% y~ equation
	a_y		= fullols(HP_cycle, [ mlag(HP_cycle,2) mean(RR_2_lags,2) ], 1);
	% pi equation
	b_pi	= fullols( inflt - mean(INFL_4_lags(:,2:4),2) , ...
				 [ ( lag(inflt) - mean(INFL_4_lags(:,2:4),2) ) lag(HP_cycle)], 1);
	% STATE 3: SET INITIAL VALUES TO BE USED (same order to match initVals of HLW)
	% Theta3 vector is: [a_y1,a_y2,a_r,b_pi,b_y, s_y~,s_pi,s_y*]
	% 									[   1,   2,  3,   4,  5,    6,   7,  8 ] 
	initVals3 = [	a_y.bhat; 
								b_pi.bhat; 
								sqrt(a_y.sig2); 
								sqrt(b_pi.sig2);
								std( diff(HP_trend) ) ];
	% Stage 3 full initivalues with g and z as well						
	initVals3_gz = [initVals3; .01; .01];
 	% state initial values as in HLW (P00 is simplified a bit, without the one period update)
	a00 = [	HP_trend([4 3 2]);      % y*(t-1) y*(t-2) y*(t-3)
					delta_HP_trend([4 3]);	% g(t-1) g(t-2) (actually, it should be g(t-2) g(t-3), but this make no diff)
					zeros(2,1)];						% z(t-1) z(t-2) 
	P00 = 0.2*eye(length(a00));
%  	P00 = 1e6*eye(length(a00));

	% set Stage 3 state initial valalues					
	other_pars_S3.a00 = a00;							
	other_pars_S3.P00 = P00;							
		
	%% ----	STAGE 3 -------------------------------------------------------------------------------------
	%		Theta3 vec.:	[a_y1,a_y2,a_r,b_pi,b_y, s_y~,s_pi,s_y*]
	%									[   1,   2,  3,   4,  5,    6,   7,  8 ] 
	% ---------------------------------------------------------------------------------------------------
	% make data input file
	data_in		= [Y X]';

	% STAGE 3: SET BOUNDS ON THE PARAMETER SPACE IF NEEDED
	% ---------------------------------------------------------------------------------------------------
	% In HLW, (b_y) b2.constraint = 0.025, (a_r) a3.constraint = -0.0025
	% ---------------------------------------------------------------------------------------------------
	S3k	= size(initVals3_gz,1);
	LB	= struct; UB	= struct; AA	= struct; 
	A_	= [	1 1; -1 1]; bb	= [.99 .99];
	% MLE SIGMA_G,SIGMA_Z
	LB.gz = -Inf*ones(S3k,1); LB.gz(5)	=	 0.0250;	% (b_y) b2.constraint >  0.0250
 	LB.gz(1) = 1.1; 
	UB.gz =  Inf*ones(S3k,1); UB.gz(3)	= -0.0025;	% (a_r) a3.constraint < -0.0025
	AA.gz	=	[A_ zeros(2,length(LB.gz)-2)];
	UB.gz(2) = -.2; % enforce second root to be negative

	% allocate structure space 
	LL3			= struct();		gradOpt3 = struct();	
	theta3	= struct();		hessOpt3 = struct();

	%% STAGE 3 MLE ESTIMATION
	% ---------------------------------------------------------------------------------------------------
	% NOTE: I use the initVals for the baseline model, and then use the baseline models estimates
	% as initial values in the estimations that follow for the paramters that are the same.
	% ---------------------------------------------------------------------------------------------------
	% LL at initival values
	LL3.initvals	= LogLike_Stage3_HLW_SSF_sigma_gz(initVals3_gz, data_in, other_pars_S3);
	
	% \SIGMA_G AND \SIGMA_Z FREELY ESTIMATED BY MLE 
	[theta3.gz, LL3.gz]	= fmincon(@LogLike_Stage3_HLW_SSF_sigma_gz,	initVals3_gz,	AA.gz, bb,[],[],LB.gz, UB.gz, [],	options_con, data_in, other_pars_S3);
	
	% take absolute vaues of standard deviations to have postive values, -> does not matter for estimation
	theta3.gz(end-4:end)	= abs(theta3.gz (end-4:end));
	% ---------------------------------------------------------------------------------------------------
	S3_initVals		= [initVals3_gz; -LL3.initvals];
% 	S3_HLW_Rfile	= [HLW_Stage3_R_File(1:11,2)];	
	S3_MLE_gz			= [theta3.gz; -LL3.gz];
	
	% PRINT RESULTS TO SCREEN / THIS MAKES A MATLABE TABLE OBJECT
	Cols2Print = 1:5;
% 	tabS3 = table(par_names3, S3_initVals, S3_HLW_Rfile, S3_MLE_gz);
	tabS3 = table(par_names3, S3_initVals, S3_MLE_gz);
	sep(120)
	fprintf('																Stage 3 Results           \n')
	sep(120)
% 	print2screen(tabS3([1:10 12:13 11],[1 3:end]), '%16.10f')
	print2screen(tabS3, '%16.10f'); sep(120)

	S3MLE = [S3MLE S3_MLE_gz];
	end

%% Now arrange and plot the series 
% clc
par_names3(6) = {'sigma_ytld'};
par_names3(8) = {'sigma_ystar'};
par_names3(11)= {'Log_Like'};
rec_pars_S3 = array2timetable(S3MLE','VariableNames',par_names3,'RowTimes', data.Time(TT) );
% head2tail(rec_pars_S3(:,[9 10]));
% save('recursive_stage3_estimates_US_diffuse.mat','rec_pars_S3')
% print2xls(rec_pars_S3,'rec_pars_S3_diffuse.xls')
load('_some_output/recursive_stage3_estimates_US_sqp.mat')

%
clc
clf
% f1 = figure(1); clf; set(f1,'WindowState','maximized','Position',[1441 641 1200 1843]);
set(groot,'defaultLineLineWidth',2.5); 
% PLOTTING DIMENSIONS FOR SUB-FIGURES ---------------------------------------------------------------
p.fs = 16;			% FONT SIZE                                       
p.df = 12;			% (INVERSE) DATE FREQUENCY, reduce to have more space between dates
p.aj = -1.225;	% SUBTITLE POSITION ADJUSTMENT                  
p.dm = @(x)([.3-(x-1)*.145 .86 .285]);% FIG DIMENSION      

plot(rec_pars_S3{:,9}); hold on; 
plot(rec_pars_S3{:,10},'--','LineWidth',3)

hold off; 
grid on; box on;
setplot(p.dm(0),p.fs,2)
setdateticks(rec_pars_S3.Time, p.df)
setyticklabels( 0:.01:.06, 2); 
ylim([-.002 .062]);
hline(0);
set(gca,'GridLineStyle',':','GridAlpha',1/3);
% set(gca,'XTickLabelRotation',0);
setoutsideTicks; add2yaxislabel
addlegend({'$\sigma_g$';'$\sigma_z$'},1,16)
% subtitle('ML estimates of $\sigma_g$ and $\sigma_z$ from full Stage 3 model without MUE',-1.11,[],1)
print2pdf('rec_estimates_US','../../graphics')

% Print tabel to screen
% tabS3 = table(par_names3, S3MLE);
% print2screen(tabS3, '%16.10f'); sep(120)
% plot(S3MLE(1,:))

% end	














disp('Done')




































% EOF end of file 
% ---------------------------------------------------------------------------------------------
% 	figure('WindowState','maximized','Position',[1441 641 1200 1843]);
% 	figure('WindowState','maximized','Position',[1 641 1200 1843]);
%		figure('WindowState','maximized','Position',[-1441 641 1200 1843]);
% ---------------------------------------------------------------------------------------------



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	