classdef dwiMRI_Session  < dynamicprops & matlab.mixin.SetGet
    %%% Written by: Aaron Schultz (aschultz@martinos.org)
    %%%             Rodrigo Perea (rpereacamargo@mgh.harvard.edu)
    
    
    properties (GetAccess=private)
        %%% For properties that cannot be accessed or changed from outside
        %%% the class.
        %test = '12345';
    end
    
    properties (Constant)
        %%% For setting properties that never change.
    end
    
    properties (Dependent)
        %%% For setting properties that get automatically updated when
        %%% other variables change
    end
    
    properties
        sessionname = '';
        rawfiles = '';
        root=''; %to be modified if specified as an argument. 
        
        %         collectiondate = '';
        %         modality = '';
        %         type = '';
        %         proj='';
        
        dosave = false;
        wasLoaded = false;
        
        rawstructural = '';
        
        %%% Change things so that these are the only version used.
        fsdir = '';
        fsubj = '';
        
        %         surfRend = [];
        %         vbmdir = '';
        %
        %         interporder = 5;
        %
        objectHome = '';
        %
        history = [];
        error_messages = '';
        
        %
        %         TR = [];
        %         nVols = [];
        %         nRuns = [];
        %
        pth = [];
        dbentry = [];
        %
        %         bb = [-78 -112 -70; 78 76 90];
        %         bb = [];
        vox = [];
        lastFN = [];
        %
        Params = [];
        %
        %         RunInfo = [];
%         Behavior = [];
    end
    
    methods 
        function obj=dwiMRI_Session()
            setDefaultParams(obj);
        end
        
        function obj = getPaths(obj,proj)
            obj.pth = MyPaths(proj);
        end
        
        function obj = getDB(obj)
            R = DataCentral(['SELECT * FROM Sessions.MRI WHERE MRI_Session_Name="' obj.sessionname '"']);
            obj.dbentry = R;
            try
                obj.collectiondate = R.MRI_SessionDate;
            end
        end
        
        function showUnpackedData(obj)
            try
                type([obj.session_location filesep 'LogFiles' filesep 'config']);
            catch
                disp('config file was not found');
            end
        end
        
        function obj = setDefaultParams(obj)
            %%%%%%%PARAMETERS SET BASED ON EACH STEP OF THE DWI PIPELINE%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For proc_dcm2nii:
            obj.Params.DCM2NII.in.first_dcmfiles = [];
            obj.Params.DCM2NII.scanlog = '';
            obj.Params.DCM2NII.seq_names='';
            obj.Params.DCM2NII.in.nvols = [];
            obj.Params.DCM2NII.in.fsl2std_matfile = '';
            obj.Params.DCM2NII.in.fsl2std_param = '';
            
            obj.Params.DCM2NII.out.location = [];
            obj.Params.DCM2NII.out.fn = [];
            obj.Params.DCM2NII.out.bvecs='';
            obj.Params.DCM2NII.out.bvals='';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.DropVols.in.dropVols = [];
            obj.Params.DropVols.in.prefix = 'dv_';
            obj.Params.DropVols.in.movefiles = '';
            obj.Params.DropVols.in.fn = [];
            obj.Params.DropVols.out.fn = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.GradNonlinCorrect.in.movefiles = '';
            obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
            obj.Params.GradNonlinCorrect.in.gradfile = '';
            obj.Params.GradNonlinCorrect.in.fn = '';
            obj.Params.GradNonlinCorrect.out.b0='';
            obj.Params.GradNonlinCorrect.in.target = '';
            obj.Params.GradNonlinCorrect.in.fslroi= [0 1];
            obj.Params.GradNonlinCorrect.out.b0='';
            obj.Params.GradNonlinCorrect.out.warpfile = [];
            obj.Params.GradNonlinCorrect.out.meannii = [];
            obj.Params.GradNonlinCorrect.out.fn = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.Bet2.in.movefiles = '';
            obj.Params.Bet2.in.fracthrsh = [];
            obj.Params.Bet2.in.prefix = 'bet2_';
            obj.Params.Bet2.in.fn = '';
            obj.Params.Bet2.out.mask = '';
            obj.Params.Bet2.out.skull = '';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.Eddy.in.movefiles = '';
            obj.Params.Eddy.in.prefix = 'eddy_';
            obj.Params.Eddy.in.fn = '';            
            obj.Params.Eddy.in.mask = '';
            obj.Params.Eddy.in.index= [];
            obj.Params.Eddy.in.acqp= [];
            obj.Params.Eddy.in.bvals='';
            obj.Params.Eddy.in.bvecs='';
            
            obj.Params.Eddy.out.fn = '';    
            obj.Params.Eddy.out.bvecs = '';    
            obj.Params.Eddy.out.fn_acqp='';
            obj.Params.Eddy.out.fn_index='';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.B0mean.in.movefiles = '';
            obj.Params.B0mean.in.fn='';
            obj.Params.B0mean.in.b0_nvols=[];
            obj.Params.B0mean.out.fn='';
            obj.Params.B0mean.out.allb0s='';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.Dtifit.in.movefiles = '';
            obj.Params.Dtifit.in.fn='';
            obj.Params.Dtifit.in.bvecs='';
            obj.Params.Dtifit.in.bvals='';
            obj.Params.Dtifit.in.mask='';
            obj.Params.Dtifit.out.prefix='';
            obj.Params.Dtifit.out.FA='';
            obj.Params.Dtifit.out.RD = '';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.GQI.in.movefiles = '';
            obj.Params.GQI.in.fn = '' ;
            obj.Params.GQI.in.bvecs = '';
            obj.Params.GQI.in.bvals = '';
            obj.Params.GQI.in.mask = '' ;
            
            obj.Params.GQI.in.method = '4'; %for gqi 
            obj.Params.GQI.in.num_fiber = '3';  %modeling 3 fiber population
            obj.Params.GQI.in.param0 = '1.25';
            %obj.Params.GQI.sh = '';
            obj.Params.GQI.out.btable = '';
            obj.Params.GQI.out.src_fn = '';
            obj.Params.GQI.out.fibs_fn = '' ;
            obj.Params.GQI.out.fibs_GFA = '';
            obj.Params.GQI.out.export =   'gfa,nqa0,nqa1';          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.AntsReg.in.movefiles = '';
            obj.Params.AntsReg.in.fn = '';
            obj.Params.AntsReg.in.ref = '';
            obj.Params.AntsReg.in.dim = '3' ;
            obj.Params.AntsReg.in.niter = '1' ;
            obj.Params.AntsReg.in.transform = 's' ;  
            obj.Params.AntsReg.in.radius = '4' ;
            obj.Params.AntsReg.in.precision = 'd' ;
            obj.Params.AntsReg.in.prefix = '';
            obj.Params.AntsReg.out.fn = '';
            obj.Params.AntsReg.out.FA = '';
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            obj.Params.Skeletonize.in.movefiles = ['..' filesep 'Skeletonize' ];
            obj.Params.Skeletonize.in.movefiles = '';
            obj.Params.Skeletonize.in.fn = '' ;
            obj.Params.Skeletonize.in.meanFA = '';
            obj.Params.Skeletonize.in.skel_dst = '';
            obj.Params.Skeletonize.in.thr = '0.3';
            obj.Params.Skeletonize.in.ref_region = '';
            obj.Params.Skeletonize.in.prefix = [ obj.sessionname '_skel' ] ;
            obj.Params.Skeletonize.out.fn = '' ;
            obj.Params.Skeletonize.out.fn_blind = '' ;
            obj.Params.Skeletonize.in.FA= '' ; 
            obj.Params.Skeletonize.out.FA = '' ;
            
%             %%%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%
%             %%ALL THESE BELOW ARE INSTANCES OF A PREVIOUS CLASS (USED FOR
%             %%CODE RECYLING...)
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.CleanData.in.movefiles = '';
%             obj.Params.CleanData.in.prefix = 'cl_';
%             obj.Params.CleanData.in.fn = [];
%             obj.Params.CleanData.in.filter = true;
%             obj.Params.CleanData.in.motion = true;
%             obj.Params.CleanData.in.physio = true;
%             obj.Params.CleanData.in.other = false;
%             obj.Params.CleanData.in.otherFiles = [];
%             obj.Params.CleanData.in.deriv = true;
%             obj.Params.CleanData.in.square = true;
%             obj.Params.CleanData.in.reduce = true;
%             obj.Params.CleanData.out.fn = [];
%             obj.Params.CleanData.out.REGS = [];
%             obj.Params.CleanData.out.indices = [];
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.SNR.in.movefiles = 'SNR_Images';
%             obj.Params.SNR.in.fn = [];
%             obj.Params.SNR.in.thresh = 115;
%             
%             obj.Params.SNR.out.mean = [];
%             obj.Params.SNR.out.sd = [];
%             obj.Params.SNR.out.snr = [];
%             obj.Params.SNR.out.report = [];
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.QA.in.movefiles = '';
%             obj.Params.QA.in.logname = 'QA_Log';
%             obj.Params.QA.in.thinga = [];
%             obj.Params.QA.in.checksnr = true;
%             obj.Params.QA.in.SNRthresh = 98; %changed from 115 - Fdu 170224
%             obj.Params.QA.in.GlobalSigThresh = 2.5;
%             obj.Params.QA.in.MeanMovementThresh = 0.5; %changed from 0.75 - Fdu 170224
%             obj.Params.QA.in.TotalMovementThresh = 5;
%             obj.Params.QA.in.TotalRotationThresh = 5;
%             obj.Params.QA.in.MoveThresh = .75;
%             obj.Params.QA.in.RotThresh = 1.5;
%             obj.Params.QA.in.BadVolThresh = 20;
%             
%             obj.Params.QA.out.SNR = [];
%             obj.Params.QA.out.meanBold = [];
%             obj.Params.QA.out.badVols = [];
%             obj.Params.QA.out.meanMV = [];
%             obj.Params.QA.out.sdMV = [];
%             obj.Params.QA.out.reasons = [];
%             obj.Params.QA.out.badruns = [];
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.Map_To_Surface.in.movefiles = '';
%             obj.Params.Map_To_Surface.in.fn = [];
%             
%             obj.Params.Map_To_Surface.out.fn = [];
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             obj.Params.SurfRend.fsubj = obj.fsubj;
%             obj.Params.SurfRend.fsdir = obj.fsdir;
%             
%             obj.Params.SurfRend.figno = 101;
%             obj.Params.SurfRend.newfig = 1;
%             obj.Params.SurfRend.input_lh = [];
%             obj.Params.SurfRend.input_rh = [];
%             obj.Params.SurfRend.overlaythresh = [0 0];
%             obj.Params.SurfRend.colorlims = [0 Inf];
%             obj.Params.SurfRend.colomap = 'jet';
%             obj.Params.SurfRend.direction =  '+';
%             obj.Params.SurfRend.reverse = 0;
%             obj.Params.SurfRend.round = 0;
%             obj.Params.SurfRend.surface = 'pi';
%             obj.Params.SurfRend.shading =  'mixed';
%             obj.Params.SurfRend.shadingrange = [-2 2];
%             obj.Params.SurfRend.Nsurfs = 2;
            
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% BEGIN Things %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        function showHistory(obj)
            disp('');
            disp('HISTORY:');
            for ii = 1:numel(obj.history)
                disp([sprintf('%0.3i',ii) '  ' regexprep(obj.history{ii}.lastRun,'\n','')]);
            end
            disp('');
        end

        function showErrors(obj)
            disp('');
            fprintf('\n\nERROR MESSAGE HISTORY:');
            if isempty(obj.error_messages)
                fprintf(' \t \t No errors found!\n\n');
            else
                for ii = 1:numel(obj.error_messages)
                    fprintf(obj.error_messages{ii,:});
                end
            end
            disp('');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% END Things %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% BEGIN Data Processing Methods %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = proc_dcm2nii(obj,fn,dest,out_filename)
            fprintf('\n\n%s\n', ['CONVERTING DCM TO NII ' num2str(numel(obj.Params.DCM2NII.in)) ' VOLUMES']);
            %%% MODULE IS COMPLETE
            wasRun = false;
            nfn = [];
            for ii=1:numel(obj.Params.DCM2NII.in)
               if nargin>2 && ~isempty(dest)
                    obj.Params.DCM2NII.out(ii).location = dest;
                end
                
                if nargin>3 && ~isempty(out_filename)
                    obj.Params.DCM2NII.out(rr).out_filename = out_filename;
                end
                
                %Shortening naming comventions:
                clear in_file out_file 
                in_file=strtrim([obj.dcm_location filesep obj.Params.DCM2NII.in(ii).first_dcmfiles]);
                out_file=strtrim([obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.out(ii).fn ]);
                outpath=obj.Params.DCM2NII.out(ii).location;
                obj.Params.DCM2NII.in.fsl2std_matfile = [outpath 'fsl2std.matfile'];
                
                obj.Params.DCM2NII.out.bvecs=[ outpath strrep(obj.Params.DCM2NII.out.fn,'.nii.gz','.voxel_space.bvecs') ];
                obj.Params.DCM2NII.out.bvals=[ outpath strrep(obj.Params.DCM2NII.out.fn,'.nii.gz','.bvals') ];
                %Create out_file directory if doesnt exist:
                if ~exist(obj.Params.DCM2NII.out(ii).location,'dir')
                    clear exec_cmd
                    exec_cmd = ['mkdir -p ' obj.Params.DCM2NII.out(ii).location ];
                    system(exec_cmd);
                end
                
                %Processing starts here:
                if exist(in_file,'file') ~= 0 %check if in_file exists
                    if exist(out_file,'file') == 0 %check if out_file exists
                        %Check whether we get the specific number of volumes:
                        if  obj.Params.DCM2NII.specific_vols == obj.Params.DCM2NII.in(ii).nvols
                            %Somehow save my cmd in here...??
                            clear exec_cmd
                            exec_cmd=['mri_convert ' in_file ' ' out_file ];
                            system(exec_cmd)
                            
                            %Reorient to std --> 
                            fprintf('\nFslreorienting to standard...')
                            if isempty(obj.Params.DCM2NII.in.fsl2std_param)
                                exec_cmd=['fslreorient2std ' out_file ' > ' obj.Params.DCM2NII.in.fsl2std_matfile ];
                                system(exec_cmd);
                            else
                                exec_cmd=['echo -e ''' obj.Params.DCM2NII.in.fsl2std_param ''' > ' obj.Params.DCM2NII.in.fsl2std_matfile ];
                                system(exec_cmd);
                            end
                            
                            exec_cmd=['fslreorient2std ' out_file ' ' out_file ];
                            system(exec_cmd);
                            
                            %Now dealing with bvecs:
                            disp('Fslreorienting the bvecs now...')
                            temp_bvec=[outpath 'temp.bvec' ] 
                            exec_cmd=[obj.rotatae_bvecs_sh ' ' ...
                                ' ' obj.Params.DCM2NII.out.bvecs ...
                                ' ' obj.Params.DCM2NII.in.fsl2std_matfile ...
                                ' ' temp_bvec  ]; 
                            system(exec_cmd);
                            
                            system(['mv ' temp_bvec ' ' obj.Params.DCM2NII.out.bvecs ]);
                            fprintf('\n....done');
                            obj.UpdateHist(obj.Params.DCM2NII,'proc_dcm2nii',out_file,wasRun);
                        else
                            error('==> obj.Params.DCM2NII.specific_vols  not equal to obj.Params.DCM2NII.in(ii).nvols ');
                        end
                    else
                        disp([ '==> out_file: ' out_file ' exists. SKIPPING...'])
                    end
                else
                    error(['Error in obj.proc_dcm2nii ==> in_file: ' in_file 'does not exist!']);
                end
            end
            %%%
            
            
            fprintf('\n');
        end
        
        function obj = proc_gradient_nonlin_correct(obj)
            wasRun = false;
            target_all = obj.Params.GradNonlinCorrect.in.target;
            
            fn{1}=target_all{1}; %this will take set4 in ADRC data (shouldn;t affect it. Need to double-check!)
            [a b c ] = fileparts(fn{1});
            outpath=obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);
            
            %first extract the first b0s (if it doesn't exist):
            obj.Params.GradNonlinCorrect.in.b0{1}=[outpath 'firstb0_' b c ];
            if exist(obj.Params.GradNonlinCorrect.in.b0{1},'file')==0
                fprintf(['\nExtracting the first b0 of: ' obj.Params.GradNonlinCorrect.in.b0{1} ]);
                exec_cmd=['fslroi ' fn{1} ' ' obj.Params.GradNonlinCorrect.in.b0{1} ...
                    ' ' num2str(obj.Params.GradNonlinCorrect.in.fslroi) ];
                system(exec_cmd);
                fprintf('...done\n');
            end
            %now creating grad-nonlinearity in first_b0s:
            obj.Params.GradNonlinCorrect.out.b0{1}=[outpath 'gnc_firstb0_' b c ];
            first_b0_infile = obj.Params.GradNonlinCorrect.in.b0{1};
            first_b0_outfile = obj.Params.GradNonlinCorrect.out.b0{1};
            gradfile = obj.Params.GradNonlinCorrect.in.gradfile;
            
            %%% Compute the grdient nonlinearity correction
            obj.Params.GradNonlinCorrect.out.warpfile{1} = strrep(first_b0_infile,'.nii','_deform_grad_rel.nii');
            if exist(obj.Params.GradNonlinCorrect.out.warpfile{1},'file')==0
                cmd=['sh ' obj.sh_gradfile ' '  first_b0_infile ' ' first_b0_outfile ' ' gradfile ' '];
                system(cmd);
                wasRun = true;
            else
                fprintf([ 'Gnc warp file:' first_b0_outfile ' exists. Skipping...\n'])
            end
            
            %%Apply the correction to the first_b0
            if exist(first_b0_outfile,'file')==0
                exec_cmd=['applywarp -i ' first_b0_infile ' -r ' first_b0_infile ...
                    ' -o ' first_b0_outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile{1} ...
                    ' --interp=spline' ];
                fprintf(['\nGNC: Applying warp to first_b0_file: '  first_b0_infile]);
                system(exec_cmd);
                fprintf(['...done\n']);
            end
            
            
            %%% Apply the correction to all the subsequent diffusion images.
            for ii=1:numel(target_all)
                dwi_infile{ii}=target_all{ii};
                [a b c ] = fileparts(dwi_infile{ii});
                dwi_outfile{ii}=[outpath 'gnc_' b c ];
                obj.Params.GradNonlinCorrect.in.fn{ii}=dwi_infile{ii};
                obj.Params.GradNonlinCorrect.out.fn{ii}=dwi_outfile{ii};
                if exist(dwi_outfile{ii},'file')==0
                    fprintf(['\nGNC: Applying warp field to all the other images: ' dwi_infile{ii}]);
                    exec_cmd = ['applywarp -i ' dwi_infile{ii} ' -r ' first_b0_infile ...
                        ' -o ' dwi_outfile{ii} ' -w ' obj.Params.GradNonlinCorrect.out.warpfile{1} ' --interp=spline'];
                    system(exec_cmd);
                    fprintf('....done\n');
                    wasRun = true;
                end
            end
            obj.Params.GradNonlinCorrect.in.fn=obj.Params.GradNonlinCorrect.in.fn';
            obj.Params.GradNonlinCorrect.out.fn=obj.Params.GradNonlinCorrect.out.fn';
            
            obj.UpdateHist(obj.Params.GradNonlinCorrect,'proc_gradient_nonlin_correct',obj.Params.GradNonlinCorrect.out.warpfile{1},wasRun);
        end
        
        function obj = proc_bet2(obj)
            for ii=1:numel(obj.Params.Bet2.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Bet2.in.fn{ii})
                    cur_fn=cell2char(obj.Params.Bet2.in.fn{ii});
                else
                    cur_fn=obj.Params.Bet2.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.Bet2.in.movefiles);
                clear outfile
                obj.Params.Bet2.out.skull{ii} = [ outpath obj.Params.Bet2.in.prefix  b c ];
                obj.Params.Bet2.out.mask{ii} = strrep(obj.Params.Bet2.out.skull{ii},'.nii.gz','_mask.nii.gz');
                if exist( obj.Params.Bet2.out.mask{ii},'file')==0
                    fprintf(['\nExtracting the skull using bet2 for : ' obj.Params.Bet2.in.fn{ii} ]);
                    exec_cmd=[ 'bet2 ' obj.Params.Bet2.in.fn{ii} ' ' obj.Params.Bet2.out.skull{ii}  ' -m -f ' num2str(obj.Params.Bet2.in.fracthrsh) ]
                    system(exec_cmd)
                    system(['mv ' obj.Params.Bet2.out.skull{ii} '_mask.nii.gz ' obj.Params.Bet2.out.mask{ii} ] ) ;
                    wasRun=true;
                    obj.UpdateHist(obj.Params.Bet2,'proc_bet2', obj.Params.Bet2.out.mask{ii},wasRun);
                else
                     fprintf(['\n proc_bet2(): ' obj.Params.Bet2.out.mask{ii} ' exists. Skipping...\n\n']) ;
                end
            end
            
        end
        
        function obj = proc_eddy(obj)
            for ii=1:numel(obj.Params.Eddy.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Bet2.in.fn{ii})
                    cur_fn=cell2char(obj.Params.Eddy.in.fn{ii});
                else
                    cur_fn=obj.Params.Eddy.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.Eddy.in.movefiles);
                clear outfile
             
                %(dependency) Attempting to create acqp file:
                obj.Params.Eddy.out.fn_acqp{ii}= [ outpath 'acqp.txt' ] ;
                if exist(obj.Params.Eddy.out.fn_acqp{ii},'file')==0
                    fprintf(['\n Creating ' obj.Params.Eddy.out.fn_acqp{ii}  ]);
                    exec_cmd=['echo " ' num2str(obj.Params.Eddy.in.acqp) ' " >> ' obj.Params.Eddy.out.fn_acqp{ii}  ];
                    system(exec_cmd);
                    fprintf(' ...done\n');
                end
                %(dependency) Attempting to create index file:
                obj.Params.Eddy.out.fn_index{ii}= [ outpath 'index.txt' ] ;
                if exist(obj.Params.Eddy.out.fn_index{ii},'file')==0
                    fprintf(['\n Creating ' obj.Params.Eddy.out.fn_index{ii}  ]);
                    exec_cmd=['echo " ' num2str(obj.Params.Eddy.in.index) ' " >> ' obj.Params.Eddy.out.fn_index{ii}  ];
                    system(exec_cmd);
                    fprintf(' ...done\n');
                end
                %Attempting to run eddy_openmp now:
                obj.Params.Eddy.out.fn{ii} = [ outpath obj.Params.Eddy.in.prefix  b c ];
                obj.Params.Eddy.out.bvecs{ii} = [ outpath obj.Params.Eddy.in.prefix strrep(b,'.nii','.eddy_rotated_bvecs') ];
                if exist( obj.Params.Eddy.out.fn{ii},'file')==0
                   try
                    fprintf(['\nApplying eddy in: ' obj.Params.Eddy.in.fn{ii} ]);
                    fprintf(['\n this will take a couple of minutes...']);
                    exec_cmd=[ 'eddy_openmp --imain=' obj.Params.Eddy.in.fn{ii} ...
                        ' --mask=' obj.Params.Eddy.in.mask{ii} ...
                        ' --index=' obj.Params.Eddy.out.fn_index{ii} ...
                        ' --acqp='  obj.Params.Eddy.out.fn_acqp{ii}  ...
                        ' --bvecs='  obj.Params.Eddy.in.bvecs{ii} ... 
                        ' --bvals=' obj.Params.Eddy.in.bvals{ii}  ...
                        ' --repol --out=' [ outpath obj.Params.Eddy.in.prefix strrep(b,'.nii','') ]  ];
                    system(exec_cmd);
                    fprintf(['...done \n']);
                    
                    wasRun=true;
                    obj.UpdateHist(obj.Params.Eddy,'proc_eddy', obj.Params.Eddy.out.fn{ii},wasRun);
                   catch
                       errormsg=['PROC_EDDY: Cannnot run eddy in: ' ...
                            obj.Params.Eddy.in.fn{ii} 'please double check parameters?'\n' ];
                       obj.UpdateErrors(errormsg);
                   end
                else
                    fprintf(['\n proc_eddy(): ' obj.Params.Eddy.out.fn{ii} ' exists. Skipping...\n\n']) ;
                end
            end
        end
        
        function obj = proc_meanb0(obj)
            for ii=1:numel(obj.Params.B0mean.in.fn)
                clear cur_fn;
                if iscell(obj.Params.B0mean.in.fn{ii})
                    cur_fn=cell2char(obj.Params.B0mean.in.fn{ii});
                else
                    cur_fn=obj.Params.B0mean.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.B0mean.in.movefiles);
                clear outfile
                %Init variable names:
                obj.Params.B0mean.out.allb0s{ii}= [ outpath 'all_b0s_' b c ] ;
                obj.Params.B0mean.out.fn{ii}= [ outpath 'b0mean_' b c ] ;
                try
                    %Attempting to create B0means: 
                    if exist( obj.Params.B0mean.out.allb0s{ii},'file')==0
                        fprintf(['\nMerging all b0s from : ' cur_fn]);
                        exec_cmd=[ 'fslroi ' cur_fn ' ' ...
                            obj.Params.B0mean.out.allb0s{ii} ' 0 ' ...
                            num2str(obj.Params.B0mean.in.b0_nvols)  ];
                        system(exec_cmd);
                        fprintf(['...done']);
                    end
                    if exist( obj.Params.B0mean.out.fn{ii},'file')==0
                        fprintf(['\n Meaning all b0s from : ' cur_fn]);
                        exec_cmd=[ 'fslmaths ' obj.Params.B0mean.out.allb0s{ii} ...
                            ' -Tmean ' obj.Params.B0mean.out.fn{ii}];
                        system(exec_cmd);
                        fprintf(['...done \n']);
                        wasRun=true;
                         obj.UpdateHist(obj.Params.B0mean,'proc_b0mean', obj.Params.B0mean.out.fn{ii},wasRun);
                    else
                        fprintf(['\n proc_b0mean(): ' obj.Params.B0mean.out.fn{ii} ' exists. Skipping...\n\n']) ;
                    end
                catch
                    errormsg=['PROC_B0MEAN: Cannnot create the following meanB0 from:'  ...
                        cur_fn 'Please check this input location!\n' ];
                    obj.UpdateErrors(errormsg);
                    
                end
            end
        end
        
        function obj = proc_dtifit(obj)
            for ii=1:numel(obj.Params.Dtifit.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Dtifit.in.fn{ii})
                    cur_fn=cell2char(obj.Params.Dtifit.in.fn{ii});
                else
                    cur_fn=obj.Params.Dtifit.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.Dtifit.in.movefiles);
                clear outfile
                %Init variable names:
                obj.Params.Dtifit.out.FA{ii} = [ outpath  obj.sessionname '_FA.nii.gz' ] ;
                obj.Params.Dtifit.out.prefix{ii} = [ outpath  obj.sessionname ];
                try
                    %Attempting to dtifit:
                    if exist( obj.Params.Dtifit.out.FA{ii},'file')==0
                        fprintf(['\nDtifit reconstruction...']);
                        exec_cmd=[ 'dtifit -k ' obj.Params.Dtifit.in.fn{ii} ... 
                            ' -o ' obj.Params.Dtifit.out.prefix{ii} ...
                            ' -m ' obj.Params.Dtifit.in.mask{ii} ...
                            ' -r ' obj.Params.Dtifit.in.bvecs{ii} ... 
                            ' -b ' obj.Params.Dtifit.in.bvals{ii} ];
                        system(exec_cmd);
                        fprintf(['...done']);
                       
                    else
                        fprintf([ obj.Params.Dtifit.out.FA{ii} ' exist.\nWas dtifit already ran? Skipping...\n'])
                    end
                catch
                    errormsg=['PROC_DTIFIT: Cannot apply dtifit:'  ...
                        'Please check dtifit input location!\n' ];
                    obj.UpdateErrors(errormsg);
                    
                end
                %Outputting RD:
                obj.Params.Dtifit.out.RD{ii} = strrep(obj.Params.Dtifit.out.FA{ii},'FA','RD');
                if exist(obj.Params.Dtifit.out.RD{ii})==0
                    
                    exec_cmd=[ ' fslmaths ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L2') ...
                        ' -add ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L3') ...
                        ' -div 2 ' obj.Params.Dtifit.out.RD{ii}  ];
                    fprintf('\nCreating RD dtifit...');
                    system(exec_cmd);
                    fprintf('...done \n');
                end
            end
        end
        
        function obj = proc_gqi(obj)
            for ii=1:numel(obj.Params.GQI.in.fn)
                try
                    clear cur_fn;
                    if iscell(obj.Params.GQI.in.fn{ii})
                        cur_fn=cell2char(obj.Params.GQI.in.fn{ii});
                    else
                        cur_fn=obj.Params.GQI.in.fn{ii};
                    end
                    [a b c ] = fileparts(cur_fn);
                    outpath=obj.getPath(a,obj.Params.GQI.in.movefiles);
                    clear outfile
                    %Init variable names:
                    obj.Params.GQI.out.btable{ii}= [ outpath  obj.sessionname '_btable.txt' ] ;
                    
                    %Attempting to create b_table:
                    if exist(obj.Params.GQI.out.btable{ii},'file')==0
                        [~, nrow ]=system(['cat ' obj.Params.GQI.in.bvecs{ii} ' | wc -l | awk  '' {print $1} '' '  ] );
                        nrow=str2num(nrow);
                        temp_bvecs{ii}=[ outpath 'temp.txt' ];
                        if nrow == 3 ; %then its in column form, change it...
                            exec_cmd=[ 'drigo_col2rows.sh ' obj.Params.GQI.in.bvecs{ii} ...
                                ' > ' temp_bvecs{ii}];
                            system(exec_cmd);
                        else
                            exec_cmd=[ 'cat ' obj.Params.GQI.in.bvecs ' >> ' temp_bvecs{ii} ];
                            system(exec_cmd);
                        end
                        exec_cmd=[' paste ' obj.Params.GQI.in.bvals{ii} ' ' ...
                            temp_bvecs{ii} ' | sed ''s/\t/ /g'' >' obj.Params.GQI.out.btable{ii}  ];
                        system(exec_cmd);
                        exec_cmd=(['rm ' temp_bvecs{ii}]);
                        system(exec_cmd);
                        %                     else
                        %                         fprintf(['\n B-table: ' obj.Params.GQI.out.btable{ii}  ' exists. Skipping creation...']);
                    end
                    
                    %Attempting to create the src.fz file:
                    obj.Params.GQI.out.src_fn{ii} = [outpath obj.sessionname '.src.gz' ];
                    if exist(obj.Params.GQI.out.src_fn{ii},'file')==0
                        fprintf(['\nSource gz file reconstruction...']);
                        exec_cmd=[ 'dsi_studio_run --action=src ' ...
                            ' --source=' obj.Params.GQI.in.fn{ii} ...
                            ' --b_table=' obj.Params.GQI.out.btable{ii} ...
                            ' --output=' obj.Params.GQI.out.src_fn{ii} ];
                        system(exec_cmd);
                        fprintf('...done');
                    else
                        fprintf(['\n The src file: ' ' exists. Skipping...\n']);
                    end
                    
                    try
                        obj.Params.GQI.out.fibs_fn{ii} = ls([outpath '*.fib.gz' ] );
                    catch
                        obj.Params.GQI.out.fibs_fn{ii} = '';
                    end
                    %Attempting to create the fib.fz file:
                    if exist(strtrim(obj.Params.GQI.out.fibs_fn{ii}),'file')==0
                        fprintf(['\nFib gz file reconstruction...']);
                        exec_cmd=[ 'dsi_studio_run --action=rec ' ...
                            ' --source=' obj.Params.GQI.out.src_fn{ii} ...
                            ' --method=' obj.Params.GQI.in.method ...
                            ' --num_fiber=' obj.Params.GQI.in.num_fiber ...
                            ' --param0=' obj.Params.GQI.in.param0 ...
                            ' --mask=' obj.Params. GQI.in.mask{ii} ];
                        system(exec_cmd);
                        fprintf('...done');
                        
                        %Assigning the fib_fn value again (if created)
                        try
                            obj.Params.GQI.out.fibs_fn{ii} = ls([outpath '*.fib.gz' ] );
                        catch
                            obj.Params.GQI.out.fibs_fn{ii} = '';
                        end
                        
                    else
                        fprintf(['\n The fib.gz file: ' obj.Params.GQI.out.fibs_fn{ii}  ' exists. Skipping...']);
                    end
                    obj.Params.GQI.out.fibs_GFA{ii} = [ strtrim(obj.Params.GQI.out.fibs_fn{ii}) '.gfa.nii.gz' ];
                    
                    %Now exporting some values (GFA,...):
                    if exist(obj.Params.GQI.out.fibs_GFA{ii},'file') == 0
                        exec_cmd=(['dsi_studio_run --action=exp ' ...
                            ' --source=' strtrim(obj.Params.GQI.out.fibs_fn{ii}) ...
                            ' --export=' obj.Params.GQI.out.export ]);
                        system(exec_cmd);
                    end
                    
                catch
                    errormsg=['PROC_GQI: Cannot complete GQI reconstruction'  ...
                        'Please check gqi input parameters (maybe sourcing dsi_studio?)!\n' ];
                    disp(errormsg)
                    obj.UpdateErrors(errormsg);
                    
                end
            end
        end
        
        function obj = proc_antsreg(obj)
            for ii=1:numel(obj.Params.AntsReg.in.fn)
                clear cur_fn;
                if iscell(obj.Params.AntsReg.in.fn{ii})
                    cur_fn=cell2char(obj.Params.AntsReg.in.fn{ii});
                else
                    cur_fn=obj.Params.AntsReg.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.AntsReg.in.movefiles);
                clear outfile
                obj.Params.AntsReg.out.fn{ii} = [ outpath obj.Params.AntsReg.in.prefix 'Warped.nii.gz' ];
                if exist(obj.Params.AntsReg.out.fn{ii},'file')==0
                    fprintf(['\nCoregistering Ants to reference... ']);
                    tic
                    exec_cmd=[ 'antsRegistrationSyN.sh ' ...
                        ' -d '  obj.Params.AntsReg.in.dim ...
                        ' -n '  obj.Params.AntsReg.in.niter ...
                        ' -t '  obj.Params.AntsReg.in.transform ...
                        ' -r '  obj.Params.AntsReg.in.radius  ...
                        ' -p '  obj.Params.AntsReg.in.precision ...
                        ' -f '  obj.Params.AntsReg.in.ref ...
                        ' -m '  obj.Params.AntsReg.in.fn{ii} ...
                        ' -o '  [ outpath obj.Params.AntsReg.in.prefix]  ];
                    system(exec_cmd)
                    time_taken=toc;
                    fprintf(['...done.\n']);
                else
                     fprintf(['\n proc_antsreg(): ' obj.Params.AntsReg.out.fn{ii} ...
                         ' exists. Skipping...\n\n']) ;
                end
                for tocomment=1:1
                    obj.Params.AntsReg.out.FA{ii} = [ outpath obj.Params.AntsReg.in.prefix 'FA.nii.gz' ];
                    if  exist(obj.Params.AntsReg.out.FA{ii},'file')==0
                        fprintf(['\n Warping dtifit metrics...']);
                        %FA:
                        exec_cmd=[ 'WarpImageMultiTransform 3 ' ...
                            ' ' obj.Params.Dtifit.out.FA{ii}  ...
                            ' ' obj.Params.AntsReg.out.FA{ii} ...
                            ' -R '  obj.Params.AntsReg.in.ref ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped','_1Warp') ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped.nii.gz','_0GenericAffine.mat')];
                        system(exec_cmd)
                        %RD:
                        exec_cmd=[ 'WarpImageMultiTransform 3 ' ...
                            ' ' obj.Params.Dtifit.out.RD{ii}  ...
                            ' ' strrep(obj.Params.AntsReg.out.FA{ii},'FA','RD') ...
                            ' -R '  obj.Params.AntsReg.in.ref ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped','_1Warp') ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped.nii.gz','_0GenericAffine.mat') ];
                        system(exec_cmd)
                        %AxD:
                        exec_cmd=[ 'WarpImageMultiTransform 3 ' ...
                            ' ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L1') ...
                            ' ' strrep(obj.Params.AntsReg.out.FA{ii},'FA','AxD') ...
                            ' -R '  obj.Params.AntsReg.in.ref ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped','_1Warp') ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped.nii.gz','_0GenericAffine.mat') ];
                        system(exec_cmd);
                        %MD:
                        exec_cmd=[ 'WarpImageMultiTransform 3 ' ...
                            ' ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','MD')  ...
                            ' ' strrep(obj.Params.AntsReg.out.FA{ii},'FA','MD') ...
                            ' -R '  obj.Params.AntsReg.in.ref ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped','_1Warp') ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped.nii.gz','_0GenericAffine.mat') ];
                        system(exec_cmd)
                        fprintf(['...done.\n']);
                    end
                end
            end 
        end
        
        function obj = proc_skeletonize(obj)
             for ii=1:numel(obj.Params.Skeletonize.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Skeletonize.in.fn{ii})
                    cur_fn=cell2char(obj.Params.Skeletonize.in.fn{ii});
                else
                    cur_fn=obj.Params.Skeletonize.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.Skeletonize.in.movefiles);
                clear outfile
                obj.Params.Skeletonize.out.fn{ii} = [ outpath obj.sessionname obj.Params.Skeletonize.in.prefix '.nii.gz' ];
                if exist(obj.Params.Skeletonize.out.fn{ii},'file')==0
                    fprintf(['\nSkeletonizing to reference... ']);
                    tic
                    exec_cmd=[ 'tbss_skeleton ' ...
                        ' -i '  obj.Params.Skeletonize.in.meanFA ...
                        ' -p '  obj.Params.Skeletonize.in.thr ...
                        ' '  obj.Params.Skeletonize.in.skel_dst ...
                        ' '  obj.Params.Skeletonize.in.ref_region  ...
                        ' '  obj.Params.Skeletonize.in.fn{ii} ...
                        ' '  obj.Params.Skeletonize.out.fn{ii} ];
                    system(exec_cmd)
                    time_taken=toc;
                    fprintf(['...done.\n']);
                else
                     fprintf(['\n proc_skeletonize(): ' obj.Params.AntsReg.out.fn{ii} ...
                         ' exists. Skipping...\n\n']) ;
                     
                end
                
                
                %NOW OTHER METRICS:
                obj.Params.Skeletonize.in.FA{ii}  = obj.Params.AntsReg.out.FA{ii} ; 
                obj.Params.Skeletonize.out.FA{ii} = [ outpath obj.sessionname obj.Params.Skeletonize.in.prefix '_FA.nii.gz' ];
                
                for tocomment=1:1;
                    if exist(obj.Params.Skeletonize.out.FA{ii},'file')==0
                        fprintf(['\nSkeletonizing  dtimetrics... ']);
                        tic
                        fprintf('\n in FA...');
                        exec_cmd=[ 'tbss_skeleton ' ...
                            ' -i '  obj.Params.Skeletonize.in.meanFA ...
                            ' -p '  obj.Params.Skeletonize.in.thr ...
                            ' '  obj.Params.Skeletonize.in.skel_dst ...
                            ' '  obj.Params.Skeletonize.in.ref_region  ...
                            ' '  obj.Params.Skeletonize.in.FA{ii} ...
                            ' '  obj.Params.Skeletonize.out.FA{ii}];
                        system(exec_cmd);
                        fprintf('...done\n');
                        %RD:
                        fprintf('\n in RD...');
                        exec_cmd=strrep(exec_cmd,'_FA.nii','_RD.nii')
                        system(exec_cmd); fprintf('...done\n');
                        %AxD:
                        fprintf('\n in AxD...');
                        exec_cmd=strrep(exec_cmd,'_RD.nii','_AxD.nii')
                        system(exec_cmd); fprintf('...done\n');
                        %MD:
                        fprintf('\n in MD...');
                        exec_cmd=strrep(exec_cmd,'_AxD.nii','_MD.nii')
                        system(exec_cmd); fprintf('...done\n');
                        toc
                        fprintf(['...done.\n']);
                    else
                        fprintf(['\n proc_skeletonize(): ' obj.Params.AntsReg.out.fn{ii} ...
                            ' exists. Skipping...\n\n']) ;
                        
                    end
                end
                %COMMENTED DUE TO UNNECESSARY USAGE (for now...):
%                 obj.getDB;                
%                 obj.Params.Skeletonize.out.fn_blind{ii} = [ outpath cell2char(obj.dbentry.SubjIDshort) obj.Params.Skeletonize.in.prefix '.nii.gz' ];
%                 if exist(obj.Params.Skeletonize.out.fn_blind{ii},'file')==0
%                     exec_cmd=[ 'cp ' obj.Params.Skeletonize.out.fn{ii} ...
%                         ' ' obj.Params.Skeletonize.out.fn_blind{ii} ];
%                     system(exec_cmd)
%                     fprintf(['...done.\n']);
%                 else
%                      fprintf(['\n proc_skeletonize(): ' obj.Params.AntsReg.out.fn{ii} ...
%                          ' exists. Skipping...\n\n']) ;
%                      
%                 end
             end
        end
            
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% END Data Processing Methods %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = protected)
        
        function [ind, assume] = CheckHist(obj,step,wasRun)
            assume = false;
            ind = [];
            
            for ii = 1:numel(obj.history)
                tmp = obj.history{ii};
                if isfield(tmp,'assume')
                    ind(ii,1) = isequaln(rmfield(tmp,'assume'),step);
                else
                    ind(ii,1) = isequaln(tmp,step);
                end
            end
            ind = find(ind==1);
            
            if isempty(ind);
                disp('Just Did');
                ind = numel(obj.history)+1;
                if ~wasRun
                    assume=true;
                end
            else
                disp('Already done');
            end
        end
        
        function obj = UpdateErrors(obj,message)
            new_msg=['\nOn: ' date() ' --> \n\t' message '\n'];
            obj.error_messages{end+1}=new_msg;
        end
        
        function obj = UpdateHist(obj,Params,cmd,checkFile,wasRun)
            if wasRun
                info = regexprep(sprintf('%-30s%s', [cmd ':'], obj.UserTime),'\n','');
            else
                %                 %%% For Linux
                [a theUser] = system(['ls -l --full-time ' checkFile ' |  awk  ''{print $3}''']);
                [a theTime] = system(['ls -l --full-time ' checkFile ' |  awk  ''{print $6" "$7}''']);
                theTime = datestr(theTime(1:end-11));
                
                %%% For Unix
                %                 [a theUser] = system(['ls -l -T ' checkFile ' |  awk  ''{print $3}''']);
                %                 [a theTime] = system(['ls -l -T ' checkFile ' |  awk  ''{print $6" "$7" "$8" "$9}''']);
                %                 theTime = datestr(theTime);
                %                 theUser = regexprep(theUser,'\n','');
                
                %info = [cmd ': Last run by ' theUser ' on ' theTime];
                info = regexprep(sprintf('%-30s%s', [cmd ':'], ['Last run by (file created) ' theUser ' on ' theTime]),'\n','');
            end
                        
            Params.lastRun = info;
            [ind assume] = obj.CheckHist(Params,wasRun);
            Params.assume = assume;
            obj.history{ind,1} = Params;
            
            if obj.dosave
                save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
            end
        end
        
        function outpath = getPath(obj,a,movefiles)
            if isempty(movefiles)
                outpath = [a filesep movefiles filesep];
            else
                if movefiles(1)==filesep || movefiles(1)=='~'
                    outpath = [movefiles filesep];
                else
                    outpath = [a filesep movefiles filesep];
                end
            end
            outpath = regexprep(outpath,[filesep filesep],filesep);
            
            %%%% Do I want to do the below?
            if ~exist(outpath,'dir');
                mkdir(outpath);
            end
            hm = pwd;
            cd(outpath)
            
            outpath = [pwd filesep];
            cd(hm);
        end
        
        function out = UserTime(obj)
            tmp = pwd;
            cd ~
            user = pwd;
            cd(tmp);
            
            ind = find(user == filesep);
            if ind(end)==numel(user);
                user = user(ind(end-1)+1:ind(end)-1);
            else
                user = user(ind(end)+1:end);
            end
            out = ['last run by ' user ' on ' datestr(clock)];
        end
        
        function obj = make_root(obj)
            if exist(obj.root,'dir')==0
                try
                    system(['mkdir -p' obj.root ]); 
                catch
                    disp([ 'Trying to create /DWIs/ in:' obj.root ... 
                        ' Maybe some permission issues?'])
                end
            end
        end
        
    end
end




