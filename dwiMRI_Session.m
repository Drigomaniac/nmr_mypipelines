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
        %         rootlocation = '';
        rawfiles = '';
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
            obj.collectiondate = R.MRI_SessionDate;
        end
        
        function viewMeanMaps(obj,opt)
            if nargin==1
                fn = dir_wfp2([obj.pth.restingdir filesep obj.sessionname filesep 'Rest*' filesep 'SNR_Images' filesep '*MEAN*.nii']);
            else
                fn = dir_wfp2([obj.pth.restingdir filesep obj.sessionname filesep 'Rest*' filesep 'SNR_Images' filesep '*' opt '*.nii']);
            end
            FIVE(fn);
        end
        
        function checkCoReg(obj)
            f1 = obj.Params.Coreg.in.target;
            f2 = obj.Params.Coreg.out.regimage;
            f3 = [fileparts(obj.Params.Coreg.in.target) '/VolumeLabels_shell.nii'];
            if exist(f3,'file')==0
                f3 = [obj.Params.Coreg.out.labels{3}(1:end-4) '_shell.nii'];
            end            
            
            if ~isempty(f1) && ~isempty(f2)
                if ~exist(f3,'file')
                    hh = FIVE({ {f1} {f2} });
                else
                    hh =FIVE({ {f1} {f2 f3} });
                end
            end
            
            set(hh(1).con(21),'Value',3); 
            set(hh(1).con(1),'Value',17);
            feval(get(hh(1).con(1),'CallBack')); shg
            
%             set(hh(1).con(21),'Value',2); 
%             set(hh(1).con(1),'Value',find(nominal(get(hh(1).con(1),'String'))=='test1'));
%             feval(get(hh(1).con(1),'CallBack')); shg
%             
%             feval(get(hh(1).con(21),'CallBack')); shg
        end
        
        function checkNorm(obj)
            f1 = obj.Params.ApplyNormNew.out.normmean;
            FIVE({ {which('avg152T1.nii')} {f1} });
        end
        
        function viewRaw4D(obj,ff)
            %fn = dir_wfp2([obj.pth.restingdir filesep obj.sessionname filesep 'Rest*' filesep obj.sessionname '*.nii']);
            if nargin==1
                fn = obj.rawfiles;
            else
                fn = ff;
            end
            
            cmd = 'fslview';
            for ii = 1:numel(fn)
                cmd = [cmd ' ' fn{ii}];
            end
            runFS([cmd ' &'],pwd);
        end
        
        function viewMotion(obj,opt)
            fn = obj.Params.Realign.out.realigpars;
            
            if nargin == 2
                mot = load(fn{opt});
            else
                mot = [];
                for ii = 1:numel(fn)
                    mot = [mot; load(fn{ii})];
                end
            end
            figure(99); clf;
            subplot(2,1,1); plot(mot(:,1:3),'o-','linewidth',2);
            if nargin ==2
                title(['Run ' sprintf('%0.2i',opt)])
            else
                title(['All Runs'])
            end
            
            xlabel('Volumes'); ylabel('translation'); legend({'x' 'y' 'z'},'location','best'); set(gca,'fontsize',14);
            subplot(2,1,2); plot(mot(:,4:6),'o-','linewidth',2);
            xlabel('Volumes'); ylabel('rotation'); legend({'pitch' 'roll' 'yaw'},'location','best'); set(gca,'fontsize',14);
        end
        
        function obj = viewSurf(obj,fn)
            
            a = obj.Params.SurfRend;            
            
            a.input_lh = fn{1};
            a.input_rh = fn{2};
            
            [h hh] = surfPlot3(a);
            
            set(hh,'FaceAlpha',.75);
            
            obj.Params.SurfRend = a;
        end
        
        function checkProc(obj,stem)
            [fn1 fn2] = dir_wfp2([obj.pth.restingdir filesep obj.sessionname filesep 'Rest*' filesep '' filesep stem]);
            disp(fn2);
        end
        
        function showCollectedData(obj)
            type([obj.pth.funcdir filesep obj.sessionname filesep 'LogFiles' filesep 'scan.log']);
        end
        
        function FS_showBMsurfs(obj)
            runFS(['tkmedit ' obj.sessionname ' brainmask.mgz -aux T1.mgz -surfs -aseg &'],obj.fsdir);
        end
        
        function FS_showWMsurfs(obj)
            runFS(['tkmedit ' obj.sessionname ' wm.mgz -aux T1.mgz -surfs &'],obj.dbentry.FreesurferLocation{1});
        end
        
        function showUnpackedData(obj)
            try
                type([obj.pth.funcdir filesep obj.sessionname filesep 'LogFiles' filesep 'config']);
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
            obj.Params.DCM2NII.in.prefix = [];
            obj.Params.DCM2NII.out.location = [];
            obj.Params.DCM2NII.out.filename = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.DropVols.in.dropVols = [];
            obj.Params.DropVols.in.prefix = 'dv_';
            obj.Params.DropVols.in.movefiles = '';
            obj.Params.DropVols.in.fn = [];
            obj.Params.DropVols.out.fn = [];
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For gradient non-linearity correction
            obj.Params.GradNonlinCorrect.in.movefiles = '../04_GradCorrect/';
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
        
        function obj = QuickICA(obj,fn)
            
            nComps = 30;
            intensityNormalization=true;
            movefiles = 'ICA/';
            % Other normalization and processing options?
            
            
            [a b c] = fileparts(fn{1});
            outpath = obj.getPath(a,movefiles);
            if ~exist(outpath,'dir'); mkdir(outpath); end
            %%%
            h = spm_vol([fn{1} ',1']);
            M = [];
            for ii = 1:numel(fn)
                %%% Establish options for these?
                % M = [M; zscore(detrend(FastRead(fn{ii})))];
                M = [M; zscore(detrend(ReadAndSmooth(fn{ii},[6 6 6],[])))];
            end
            %%%
            tmp = std(M);
            % mask_ind = find(~isnan(tmp) & tmp~=0);

            msk = openIMG(obj.Params.Coreg.out.labels{3});
            mask_ind = find(msk(:)'>0 & ~isnan(tmp) & tmp~=0);

%             data = demean(M(:,mask_ind)');
%             % Set comps based on variance?  How to identify noise/nuisance components?
            [b, dewhiteM1, Lambda1, V1, whiteM1] = icatb_calculate_pca( (M(:,mask_ind)'), nComps*2, 'remove_mean', 0);
            
            ICA_Options = {
                'block'         [153] ...
                'stop'          [1.0000e-06] ...
                'weights'       [0] ...
                'lrate'         [0.005]  ...
                'maxsteps'      [512] ...
                'anneal'        [0.9000]   ...
                'annealdeg'     [60]  ...
                'momentum'      [0]  ...
                'extended'      [0]  ...
                'ncomps'        [nComps]  ...
                'posact'        'off' ...
                'sphering'      'on'  ...
                'bias'          'on'  ...
                'verbose'       'on'};
            
            [icaAlgo, W, AA, icasig] = icatb_icaAlgorithm('Infomax', b', ICA_Options);
            %%%
            ic = W*b';
            tc = dewhiteM1*pinv(W);
            %%%
            for compNum = 1:nComps
                
                skew = icatb_skewness(ic(compNum,:));
                ic(compNum,:)=ic(compNum,:)*sign(skew);
                tc(:,compNum) = tc(:,compNum)*sign(skew);
            end
            save([outpath '/Regs.mat'],'tc');
            
            % XX = zscore(tc);
            % save([outpath '/Regs.mat'],'XX');
            %%%
            % b = pinv(XX)*data';
            
            hh = h(1);
            hh.dt = [16 1];
            vol = zeros(hh.dim);
            for ii = 1:size(ic,1)
                
                %vol(mask_ind) = b(ii,:);
                vol(mask_ind) = ic(ii,:);
                hh.fname = [outpath 'comp_' sprintf('%0.3i',ii) '.nii'];               
                spm_write_vol(hh,vol);
                
                
                % bb = pinv(zscore(XX(:,ii)))*data';
                % vol(mask_ind) = bb;
                % hh.fname = [outpath '/comp_' sprintf('%0.3i',ii) '_b.nii'];
                % spm_write_vol(hh(1),vol);
                %%% Write ICs?
                %     vol(:) = 0;
                %     vol(mask_ind) = ic(ii,:);
                %     spm_write_vol(hh(1),vol);
            end   
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% END Things %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% BEGIN Data Processing Methods %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = proc_dcm2nii(obj,fn,dest,out_filename) %not sure why fn??
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
                out_file=strtrim([obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.out(ii).filename ]);
                
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%AARON METHODS RECYCLED FROM fMRI_Session BELOW THIS LINE %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function obj = proc_slice_time(obj,fn)
            fprintf('\n\n%s\n', 'PERFORMING SLICE TIME CORRECTION:');
            %%% MODULE IS COMPLETE
            obj.Params.SliceTime.in.timing = [(obj.TR/numel(obj.Params.SliceTime.in.sliceorder)) (obj.TR/numel(obj.Params.SliceTime.in.sliceorder))-((obj.TR/numel(obj.Params.SliceTime.in.sliceorder))/numel(obj.Params.SliceTime.in.sliceorder))];
            wasRun = false;
            
            nfn = [];
            for ii = 1:length(fn)
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.Params.SliceTime.in.movefiles);
                if exist(outpath,'dir')==0
                    mkdir(outpath);
                end
                
                nfn{ii,1} = regexprep([outpath obj.Params.SliceTime.in.prefix b c],[filesep filesep],filesep);
                
                if exist(nfn{ii},'file')>0
                    disp('slice time correction is complete');
                    continue
                end
                
                wasRun = true;
                
                spm_slice_timing(fn{ii}, obj.Params.SliceTime.in.sliceorder, obj.Params.SliceTime.in.refslice, obj.Params.SliceTime.in.timing, obj.Params.SliceTime.in.prefix);
                try movefile([a filesep  obj.Params.SliceTime.in.prefix b c],nfn{ii}); end
                disp('slice time correction is complete');
            end
            
            
            obj.Params.SliceTime.in.fn = fn;
            obj.Params.SliceTime.out.fn = nfn;
            obj.lastFN = nfn;
            
            obj.UpdateHist(obj.Params.SliceTime,'proc_slice_time',nfn{end},wasRun);
            
            fprintf('\n');
        end
        
        function obj = proc_realigfn(obj,fn)
            fprintf('\n\n%s\n', 'REALIGNING IMAGES:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            rps = []; check = 0;
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.Params.Realign.in.movefiles);
                rps{jj,1} = [outpath 'rp_' b '.txt'];
                check = check+(exist(rps{jj},'file')>0);
            end
            
            
            
            if check==numel(fn)
                disp('realignment is complete');
            else
                wasRun = true;
                spm_realign(fn,obj.Params.Realign.in.pars);
                for jj = 1:numel(fn)
                    [a b c] = fileparts(fn{jj});
                    outpath = obj.getPath(a,obj.Params.Realign.in.movefiles);
                    if exist(outpath,'dir')==0; mkdir(outpath); end
                    try movefile([a filesep 'rp_' b '.txt'],[outpath 'rp_' b '.txt']); catch; end
                    rps{jj,1} = [outpath 'rp_' b '.txt'];
                end
                disp('realignment is complete');
            end
            
            obj.Params.Realign.out.realigpars=rps;
            obj.lastFN = fn;
            
            obj.UpdateHist(obj.Params.Realign,'proc_realign',rps{1},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_reslice(obj,fn)
            fprintf('\n\n%s\n', 'RESLICE IMAGES:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            
            [a b c] = fileparts(fn{1});
            outpath = obj.getPath(a,obj.Params.Reslice.in.movefiles);
            
            for jj = 1:numel(fn)
                if obj.Params.Reslice.in.pars.which>0
                    [a b c] = fileparts(fn{jj});
                    outpath = obj.getPath(a,obj.Params.Reslice.in.movefiles);
                    nfn{jj,1} = [outpath obj.Params.Reslice.in.pars.prefix b c];
                else
                    nfn{jj,1} = fn{jj};
                end
            end
            
            [a b c] = fileparts(fn{1});
            outpath = obj.getPath(a,obj.Params.Reslice.in.movefiles);
            
            if exist([outpath 'mean' b c],'file')>0
                %disp(['found ' outpath 'mean' b c]);
                obj.Params.Reslice.out.meanimage = [outpath 'mean' b c];
                disp('reslicing is complete');
            else
                wasRun = true;
                spm_reslice(fn,obj.Params.Reslice.in.pars);
                try movefile([a filesep 'mean' b c],[outpath 'mean' b c]); end
                obj.Params.Reslice.out.meanimage = [outpath 'mean' b c];
                for ii = 1:numel(fn)
                    [a b c] = fileparts(fn{ii});
                    outpath = obj.getPath(a,obj.Params.Reslice.in.movefiles);
                    try movefile([a filesep  obj.Params.Reslice.in.pars.prefix b c],[outpath obj.Params.Reslice.in.pars.prefix b c]); end
                end
                disp('reslicing is complete');
            end
            
            obj.Params.Reslice.in.fn = fn;
            obj.Params.Reslice.out.fn = nfn;
            obj.lastFN = nfn;
            
            obj.UpdateHist(obj.Params.Reslice,'proc_reslice',nfn{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_normalize_old(obj)
            fprintf('\n\n%s\n', 'NORMALIZING IMAGES:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            
            [a b c] = fileparts(obj.P_Reslice.meannii);
            outpath = obj.getPath(a,obj.P_NormOld.movefiles);
            if exist(outpath,'dir')==0
                mkdir(outpath);
            end
            
            if isempty(obj.P_NormOld.source)
                obj.P_NormOld.source=obj.P_Reslice.meannii;
            end
            
            defs = regexprep([outpath 'deformations_' b '.mat'],[filesep filesep],filesep);
            if exist(defs,'file')>0
                disp('normalization is complete');
                obj.P_NormOld.regfile = defs;
            else
                wasRun = true;
                spm_normalise(obj.P_NormOld.template,obj.P_NormOld.source,defs,'','',obj.P_NormOld);
                obj.P_NormOld.regfile = defs;
                disp('normalization is complete');
            end
            
            obj.UpdateHist(obj.P_NormOld,'proc_normalize_old',obj.P_NormOld.regfile ,wasRun);
            fprintf('\n');
        end
        
        function obj = proc_applynorm_old(obj,fn)
            fprintf('\n\n%s\n', 'APPLY OLD NORM:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            
            nfn = [];
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.P_NormOld.movefiles);
                nfn{jj,1} = regexprep([outpath obj.P_ApplyNormOld.prefix b c],[filesep filesep],filesep);
                
                if exist(nfn{jj},'file')>0
                    disp('normalization has been applied');
                else
                    wasRun = true;
                    spm_write_sn(fn{jj}, obj.P_NormOld.regfile, obj.P_ApplyNormOld);
                    try movefile([a filesep obj.P_ApplyNormOld.prefix b c],[outpath obj.P_ApplyNormOld.prefix b c]); end
                    disp('normalization has been applied');
                end
            end

            [a b c] = fileparts(obj.P_NormOld.source);
            outpath = obj.getPath(a,obj.P_NormOld.movefiles);
            if exist([outpath obj.P_ApplyNormOld.prefix b c],'file')==0
                wasRun = true;
                spm_write_sn(obj.P_NormOld.source, obj.P_NormOld.regfile, obj.P_ApplyNormOld);
                try movefile([a filesep obj.P_ApplyNormOld.prefix b c],[outpath obj.P_ApplyNormOld.prefix b c]); end
            end

            obj.P_ApplyNormOld.in = fn;
            obj.P_ApplyNormOld.out = nfn;
            obj.lastFN = nfn;
            
            obj.UpdateHist(obj.P_ApplyNormOld,'proc_applynorm_old',nfn{end} ,wasRun);
            fprintf('\n');
        end
        
        function obj = proc_smooth(obj,fn)
            fprintf('\n\n%s\n', 'SMOOTHING IMAGES:');
            %%% MODULE IS COMPLETE
            wasRun = false;

            [a b c] = fileparts(fn{1});
            outpath = obj.getPath(a,obj.Params.Smooth.in.movefiles);
            
            nfn = [];
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.Params.Smooth.in.movefiles);
                nfn{jj,1} = [outpath obj.Params.Smooth.in.prefix b c];
            end
            
            h = spm_vol([fn{1} ',1']);
            if h.dt(1)<16
                res = 16;
            else
                res=0;
            end
            
            for ii = 1:length(fn)
                if exist(nfn{ii},'file')==0
                    wasRun = true;
                    [a b c] = fileparts(fn{ii});
                    spm_smooth(fn{ii},nfn{ii},obj.Params.Smooth.in.kernel,res);
                end
                disp('spatial smoothing is complete.');
            end
            
            obj.Params.Smooth.in.fn = fn;
            obj.Params.Smooth.out.fn = nfn;
            obj.lastFN = nfn;
            
            obj.UpdateHist(obj.Params.Smooth,'proc_smooth',nfn{end} ,wasRun);
            fprintf('\n');
        end
        
        function obj = proc_normalize_new(obj,fn)
            fprintf('\n%s\n', 'NORMALIZING IMAGES USING SPM12 METHOD:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            if nargin == 1
                fn = obj.Params.NormNew.in.source;
            else
                obj.Params.NormNew.in.source=fn;
            end
            
            if iscell(fn);
                fn = char(fn);
            end
            
            
            [a b c] = fileparts(fn(1,:));
            outpath = obj.getPath(a,obj.Params.NormNew.in.movefiles);
            
            if exist(outpath,'dir')==0
                mkdir(outpath);
            end
            
            reg = regexprep([outpath 'y_' b c],[filesep filesep],filesep);
            if exist(reg,'file')>0
                disp('Deformation fields have already been computed');
                obj.Params.NormNew.out.regfile = reg;
                obj.Params.NormNew.out.iregfile = regexprep([outpath 'iy_' b c],[filesep filesep],filesep);
                
                obj.Params.NormNew.out.estTPM{1,1} = regexprep([outpath 'c1' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{2,1} = regexprep([outpath 'c2' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{3,1} = regexprep([outpath 'c3' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{4,1} = regexprep([outpath 'c4' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{5,1} = regexprep([outpath 'c5' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{6,1} = regexprep([outpath 'c6' b c],[filesep filesep],filesep);
            else
                wasRun = true;
                obj.Params.NormNew.in.nn.channel.vols = {fn};
                spm_preproc_run(obj.Params.NormNew.in.nn);
                
                
                try
                    movefile([a filesep 'c1*.nii'],outpath);
                    movefile([a filesep 'c2*.nii'],outpath);
                    movefile([a filesep 'c3*.nii'],outpath);
                    movefile([a filesep 'c4*.nii'],outpath);
                    movefile([a filesep 'c5*.nii'],outpath);
                    movefile([a filesep 'c6*.nii'],outpath);
                    movefile([a filesep 'iy*.nii'],outpath);
                    movefile([a filesep 'y*.nii'],outpath);
                    movefile([a filesep '*seg8.mat'],outpath);
                end
                
                obj.Params.NormNew.out.regfile = reg;
                obj.Params.NormNew.out.iregfile = regexprep([outpath 'iy_' b c],[filesep filesep],filesep);
                
                obj.Params.NormNew.out.estTPM{1,1} = regexprep([outpath 'c1' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{2,1} = regexprep([outpath 'c2' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{3,1} = regexprep([outpath 'c3' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{4,1} = regexprep([outpath 'c4' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{5,1} = regexprep([outpath 'c5' b c],[filesep filesep],filesep);
                obj.Params.NormNew.out.estTPM{6,1} = regexprep([outpath 'c6' b c],[filesep filesep],filesep);
            end
            
            obj.UpdateHist(obj.Params.NormNew,'proc_normalize_new',reg,wasRun);
            
        end
        
        function obj = proc_applynorm_new(obj,fn,regfile)
            fprintf('\n%s\n', 'APPLYING SPM12 STYLE NORMALIZATION:');
            %%% MODULE IS COMPLETE
            wasRun = false;

            if nargin < 2
                fn = obj.Params.ApplyNormNew.in.fn;
            else
                obj.Params.ApplyNormNew.in.fn=fn;
            end
            
            if nargin < 3
                regfile = obj.Params.ApplyNormNew.in.regfile;
            else
                obj.Params.ApplyNormNew.in.regfile=regfile;
            end
            
            if ischar(fn)
                fn = cellstr(fn);
            end            
                
            [a b c] = fileparts(fn{1});
            outpath = obj.getPath(a,obj.Params.ApplyNormNew.in.movefiles);
            ff = [outpath obj.Params.ApplyNormNew.in.prefix b c];
            
            nfn = [];
            for ii = 1:numel(fn);
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.Params.ApplyNormNew.in.movefiles);
                if exist(outpath,'dir')==0; mkdir(outpath); end
                
                ff = regexprep([outpath obj.Params.ApplyNormNew.in.prefix b c],[filesep filesep],filesep);
                nfn{end+1,1} = ff;
            end
                        
            check = 0;
            for ii = 1:numel(nfn)
                if exist(nfn{ii},'file')>0
                    check = check+1;
                end
            end
            
            if check ~= (numel(fn))
                wasRun = true;
                defs = obj.Params.ApplyNormNew.in.pars;
                defs.comp{1}.def = {regfile};
                defs.out{1}.pull.fnames  = fn(:)';
                
                spm_deformations(defs);               
                
                nfn = [];
                for ii = 1:numel(fn);
                    [a b c] = fileparts(fn{ii});
                    outpath = obj.getPath(a,obj.Params.ApplyNormNew.in.movefiles);
                    ff = [outpath obj.Params.ApplyNormNew.in.prefix b c];
                    movefile([a filesep 'w' b c], ff);
                    nfn{end+1,1} = ff;
                end
                
            end
            
            %[a b c] = fileparts(obj.Params.Reslice.out.meanimage);
            %outpath = obj.getPath(a,obj.Params.ApplyNormNew.in.movefiles);
            %if exist(outpath,'dir')==0; mkdir(outpath); end
%             if ~exist([outpath obj.Params.ApplyNormNew.in.prefix b c])
%                 defs = obj.Params.ApplyNormNew.in.pars;
%                 defs.comp{1}.def = {regfile};
%                 defs.out{1}.pull.fnames  = {obj.Params.Reslice.out.meanimage};
%                 
%                 spm_deformations(defs);  
%                 
%                 ff = [outpath obj.Params.ApplyNormNew.in.prefix b c];
%                 movefile([a filesep 'w' b c], ff);
%                 
%                 obj.Params.ApplyNormNew.out.normmean = ff;
%             else
%                 ff = [outpath obj.Params.ApplyNormNew.in.prefix b c];
%                 obj.Params.ApplyNormNew.out.normmean = ff;
%             end
            
            disp('Applying Normalization is complete');
            obj.Params.ApplyNormNew.out.fn = nfn;
            
            obj.UpdateHist(obj.Params.ApplyNormNew,'proc_applynorm_new',nfn{1},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_apply_reservsenorm_new(obj,fn,regfile,targ)
            fprintf('\n%s\n', 'APPLYING SPM12 STYLE REVERSE NORMALIZATION:');
            wasRun = false;
            
            if nargin == 1
                fn = obj.Params.ApplyReverseNormNew.in.fn;
            else
                obj.Params.ApplyReverseNormNew.in.fn=fn;
            end
            
            if nargin <3
                regfile = obj.Params.ApplyReverseNormNew.in.regfile;
            else
                obj.Params.ApplyReverseNormNew.in.regfile=regfile;
            end
            
             if nargin <4
                targ = obj.Params.ApplyReverseNormNew.in.targ;
            else
                obj.Params.ApplyReverseNormNew.in.targ=targ;
            end
                        
%             if isempty(obj.Params.ApplyReverseNormNew.in.regfile)
%                 regfile = ''
%             end

            if ischar(fn)
                fn = cellstr(fn);
            end
            
            [a b c] = fileparts(obj.Params.ApplyReverseNormNew.in.regfile);
            outpath = obj.getPath(a,obj.Params.ApplyReverseNormNew.in.movefiles);
            if exist(outpath,'dir')==0
                mkdir(outpath);
            end

            nfn = [];
            for ii = 1:numel(fn);
                [a1 b1 c1] = fileparts(fn{ii});
                ff = regexprep([outpath obj.Params.ApplyReverseNormNew.in.prefix b1 c1],[filesep filesep],filesep);
                nfn{end+1,1} = ff;
            end
            
            check = 0;
            for ii = 1:numel(nfn)
                if exist(nfn{ii},'file')>0
                    check = check+1;
                end
            end            
            
            if check ~= (numel(fn))    
                wasRun = true;
                h = spm_vol(targ);
                x = spm_imatrix(h.mat);
            
                defs = obj.Params.ApplyNormNew.in.pars;
                defs.comp{1}.def = {regfile};
                defs.out{1}.pull.fnames  = fn(:)';
                defs.comp{2}.idbbvox.vox = abs(x(7:9));
                defs.comp{2}.idbbvox.bb = world_bb(h);
                
                spm_deformations(defs);
                
                nfn = [];
                for ii = 1:numel(fn);
                    [a1 b1 c1] = fileparts(fn{ii});
                    ff = [outpath obj.Params.ApplyReverseNormNew.in.prefix b1 c1];
                    movefile([a1 filesep 'w' b1 c1], ff);
                    nfn{end+1,1} = ff;
                end
            end
            disp('Applying Normalization is complete');
            
            obj.Params.ApplyReverseNormNew.in.fn=fn;
            obj.Params.ApplyReverseNormNew.out.fn=nfn;
            
            obj.UpdateHist(obj.Params.ApplyReverseNormNew,'proc_apply_reservsenorm_new',nfn{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_snr(obj,fn)
            fprintf('\n\n%s\n', 'ESTIMATING TEMPORAL SNR:');
            
            wasRun = false;

            for ii = 1:numel(fn)
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.Params.SNR.in.movefiles);
                 
                if exist([outpath b '_MEANimage' c],'file')==0
                    wasRun = true;
                    Data_QC(fn{ii});
                end
                
                obj.Params.SNR.out.mean{ii,1} =   [outpath  b '_MEANimage' c];
                obj.Params.SNR.out.sd{ii,1} =     [outpath b '_SDimage' c];
                obj.Params.SNR.out.snr{ii,1} =    [outpath b '_SNRimage' c];
                obj.Params.SNR.out.report{ii,1} = [outpath 'Global_SNR_Report_' b '.txt' ];
            end
            obj.Params.SNR.in.fn = fn;
            
            %disp('snr maps and snr values have been computed');
            obj.UpdateHist(obj.Params.SNR,'proc_snr',obj.Params.SNR.out.report{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_screen(obj,fn)
            disp('This section still needs to be setup');
            return
            fprintf('\n%s\n', 'SCREEN DATA FOR QUALITY:');
            obj.P_ScreenData.logfile;
            
            for ii = 1:numel(fn)
                %%
                MM = FastRead(fn{ii});
                indd = find(nanmean(MM)>(nanmean(MM(:))/8));
                GM = nanmean(MM(:,indd),2);
                
                GM1 = diff(GM);
                ind1 = find(abs(GM1)>(std(GM1)*obj.P_ScreenData.globsigthresh ))+1;
                keyboard
                %                 figure(10); clf;
                %                 plot(GM,'o-');
                %                 hold on;
                %                 plot(ind1,GM(ind1),'r*');
                %%
                
                dat = load(obj.P_Realign.realigpars{ii});
                mv = (sqrt(sum(diff(dat(:,1:3)).^2,2)));
                obj.P_ScreenData.meanmovement = mean(mv);
                obj.P_ScreenData.sdmovement = std(mv);
                
                
            end
            
            
            
            
            
            tm = DistMat(dat(:,1:3),0);
            tm = max(tm(:));
            tr = DistMat(dat(:,4:6),0);
            tr = max(tr(:));
            
            
            S = {''; ['Total Movement = ' num2str(round(1000*(tm))/1000) ' milimeters']; ['Total Rotation = ' num2str(round(1000*(tr))/1000) ' degrees']};
            if tm > P.bv.TotalMovementThresh
                disp('Junking Rest Run due to too much overall movement.');
                S{end+1} = 'Junking Rest Run due to too much overall movement.';
                save NoGo_Bad_Movement.mat nnn;
            end
            if tr > P.bv.TotalRotationThresh
                disp(['Junking Rest Run due to too much overall rotation.']);
                S{end+1} = ['Junking Rest Run due to too much overall rotation.'];
                save NoGo_Bad_Rotation.mat nnn;
            end
            S{end+1} = '';
            disp(S{2});
            disp(S{3});
            WriteDataToText(S, nnn, 'a', '\t');
            
            tmp = diff(dat).^2;
            pos = sqrt(sum(tmp(:,1:3),2));
            ori = sqrt(sum(tmp(:,4:6),2)) * (360/(2*pi));
            ind2 = find(pos>P.bv.MoveThresh)+1;
            ind3 = find(ori>P.bv.RotThresh)+1;
            
            ind = unique([ind1'; ind2; ind3]);
            R =  zeros(size(MM,4),numel(ind));
            for qq = 1:length(ind);
                R(ind(qq),qq) = 1;
            end
            
            
            S = [];
            S{1} = ['Negating ' num2str(length(ind)) ' volumes; With local indices of: '  num2str(ind')];
            WriteDataToText(S, nnn, 'a', '\t');
            
            if numel(ind)>P.bv.BadVolThresh
                disp(['Junking Rest Run due to too many bad volumes.']);
                S = [];
                S.a1 = ['Junking Rest Run due to too many bad volumes.'];
                WriteDataToText(S, nnn, 'a', '\t');
                save NoGo_Too_Many_Bad_Vols.mat nnn;
            end
            
            S = [];
            WriteDataToText({' '}, nnn, 'a', '\t');
            type(nnn);
            
            save([P.bv.BadVolRegsName(1:end-4) fn(end-4) '.mat'], 'R')
            
            %%% Finish up.
            P1.P.bv = P.bv;
            P1.P.bv.UserTime = UserTime;
            P = P1.P;
            save([pname '.mat'], 'P');
        end
        
        function obj = proc_t1_spm(obj)
            fprintf('\n\n%s\n', 'PROCESS T1 WITH SPM12:');
            wasRun = false;
            ed
            if isempty(obj.Params.spmT1_Proc.in.outdir) || isempty(obj.Params.spmT1_Proc.in.t1)
                if exist([obj.fsdir filesep obj.fsubj '/mri/spm12/'],'dir')==0
                    mkdir([obj.fsdir filesep obj.fsubj '/mri/spm12/'])
                end
                if exist([obj.fsdir filesep obj.fsubj '/mri/spm12/orig.nii'],'file')==0
                    runFS(['mri_convert ' obj.fsdir filesep obj.fsubj '/mri/orig.mgz ' obj.fsdir filesep obj.fsubj '/mri/spm12/orig.nii'],obj.fsdir);
                end
                obj.Params.spmT1_Proc.in.outdir = [obj.fsdir filesep obj.fsubj '/mri/spm12/'];
                obj.Params.spmT1_Proc.in.t1 = [obj.fsdir filesep obj.fsubj '/mri/spm12/orig.nii'];
            end
            
            root = obj.Params.spmT1_Proc.in.outdir;
            root = regexprep(root,'//','/');
            
            t1 =  obj.Params.spmT1_Proc.in.t1;
            [a b c] = fileparts(t1);

            %%%
            if exist([root filesep 'Affine.mat'],'file')==0
                %%% could put something in here to redo subseuqent steps if
                %%% this step is not complete.
                wasRun=true;
                
                in = obj.Params.spmT1_Proc.in.pars;
                in.P = t1;
                
                in.tpm = spm_load_priors8(obj.Params.spmT1_Proc.in.tpm);
                
                [Affine,h] = spm_maff8(in.P,in.samp,in.fwhm1,in.tpm,in.M,in.regtype);
                in.M = Affine;
                [Affine,h] = spm_maff8(in.P,in.samp,in.fwhm2,in.tpm,in.M,in.regtype);
                save([a filesep 'Affine.mat'],'Affine');
            else
                load([a filesep 'Affine.mat']);
            end
            fprintf('%s\n', 'Affine Registration is complete:');
            
            
            if exist([a filesep 'Norm_Seg8.mat'],'file')==0
                wasRun=true;
                
                NormPars = obj.Params.spmT1_Proc.in.NormPars;
                NormPars.image = spm_vol(t1);
                NormPars.Affine = Affine;
                NormPars.tpm = spm_load_priors8(obj.Params.spmT1_Proc.in.tpm);
                
                results = spm_preproc8(NormPars);
                save([a filesep 'Norm_Seg8.mat'],'results');
            else
                load([a filesep 'Norm_Seg8.mat']);
            end
            fprintf('%s\n', 'Normalization computation is complete:');
            
            c = regexprep(c,',1','');
            
            if exist([a filesep 'y_' b c],'file')==0
                %%% I forget what some of these static options do.  Will leave
                %%% them fixed as is for now.
                wasRun=true;
                [cls,M1] = spm_preproc_write8(results,obj.Params.spmT1_Proc.in.rflags.writeopts,[1 1],[1 1],0,1,obj.Params.spmT1_Proc.in.rflags.bb,obj.Params.spmT1_Proc.in.rflags.vox);
            end
            fprintf('%s\n', 'Normalization files have been written out:');

            if exist([a filesep 'w' b c],'file')==0
                wasRun=true;
                defs = obj.Params.spmT1_Proc.in.defs;
                defs.comp{1}.def = {[a filesep 'y_' b c]};
                defs.out{1}.pull.fnames = {t1};
                spm_deformations(defs);
            end
            obj.Params.spmT1_Proc.out.normT1 = [a filesep 'w' b c];
            fprintf('%s\n', 'Normalization has been applied to the T1:');
            
            obj.Params.spmT1_Proc.out.estTPM{1,1} = [a filesep 'c1' b c];
            obj.Params.spmT1_Proc.out.estTPM{2,1} = [a filesep 'c2' b c];
            obj.Params.spmT1_Proc.out.estTPM{3,1} = [a filesep 'c3' b c];
            obj.Params.spmT1_Proc.out.estTPM{4,1} = [a filesep 'c4' b c];
            obj.Params.spmT1_Proc.out.estTPM{5,1} = [a filesep 'c5' b c];
            obj.Params.spmT1_Proc.out.estTPM{6,1} = [a filesep 'c6' b c];
            
            obj.Params.spmT1_Proc.out.regfile = [a filesep 'y_' b c];
            obj.Params.spmT1_Proc.out.iregfile = [a filesep 'iy_' b c];
            
            obj.UpdateHist(obj.Params.spmT1_Proc,'proc_t1_spm',[a filesep 'y_' b c] ,wasRun);
            fprintf('\n');
        end
        
        function obj = proc_t1_fs(obj,fn)
            disp('Leave this to the T1 object scripts');
            return
            
%             if numel(fn)==1
%                 cmd = ['recon-all -all -nuintensitycor-3T -subjid ' obj.sessionname ' -i ' fn{1} ' &'];
%             else
%                 cmd = ['recon-all -all -nuintensitycor-3T -subjid ' obj.sessionname ' -i ' fn{1} ' -i ' fn{2} ' &'];
%             end
%             
%             if exist([obj.fsdir filesep obj.sessionname filesep 'stats' filesep 'lh.aparc.stats'],'file')>0
%                 disp('FS recon is complete');
%                 return
%             elseif exist([obj.fsdir filesep obj.sessionname],'dir')>0
%                 disp('FS folder exists, but has not completed running');
%                 return
%             else
%                 runFS(cmd,obj.fsdir);
%             end
        end
     
        function obj = proc_coreg(obj,source,target)
            fprintf('\n%s\n', 'PERFORMING COREGISTRATION');
            wasRun = false;

            if nargin<2
                source = obj.Params.Coreg.in.source;
            else
                obj.Params.Coreg.in.source = source;
            end
            
            if iscell(source)
                source = char(source);
            end
            
            if nargin<3 
                target = obj.Params.Coreg.in.target;
            else
                obj.Params.Coreg.in.target = target;
            end
            
            if iscell(target)
                target = char(target);
            end
            
            [a b c] = fileparts(source);
            outpath1 = obj.getPath(a,obj.Params.Coreg.in.movefiles);
            outpath2 = obj.getPath(a,obj.Params.Coreg.in.movecrfiles);
            
            switch lower(obj.Params.Coreg.in.style)
                case 'spm'
                    disp('USING SPM:');
                    if exist([outpath2 'CoReg.mat'],'file')==0
                        wasRun=true;
                        
                        VG = spm_vol(target);
                        VF = spm_vol(source);
                        
                        x = spm_coreg(VG,VF,obj.Params.Coreg.in.spm);
                        M  = spm_matrix(x);
                        save([outpath2 'CoReg.mat'],'x','M');
                    else
                        load([outpath2 'CoReg.mat']);
                    end
                    obj.Params.Coreg.out.regfile = [outpath2 'CoReg.mat'];
                    disp('Coregistration Estimation has been performed:');
            
                    [a1 b1 c1] = fileparts(source);
                    new = [outpath1  obj.Params.Coreg.in.prefix b1 c1];
                    if exist(new,'file')==0
                        wasRun=true;
                        copyfile(source,new);
                        MM = spm_get_space(new);
                        spm_get_space(new, M\MM);
                        if obj.Params.Coreg.in.reslice
                            h1 = spm_vol(obj.Params.Coreg.in.target);
                            h2 = spm_vol(new);
                            m = resizeVol2(h2,h1,obj.Params.Coreg.in.resampopt);
                            [a b c] = fileparts(new);
                            h1.fname = [a filesep 'rs_' b c];
                            h1.dt = h2.dt;
                            spm_write_vol(h1,m);
                        end
                    end
                    obj.Params.Coreg.out.regimage = new;
                    disp('Coregistration has been applied to the mean image:');
                    
                case 'bbreg'
                    disp('    USING BBREGISTER:');
                    [a1 b1 c1] = fileparts(source);
                    
                    reg = [outpath2 'reg_' b1 '.dat'];
                    if exist(reg,'file')==0
                        wasRun = true;
                        cmd = ['bbregister --s  '  obj.fsubj ' --mov ' source ' --init-' obj.Params.FS_Params.bbreg ' --' obj.Params.FS_Params.conweight ' --reg ' reg];
                        runFS(cmd,obj.fsdir);
                    end
                    disp('Coregistration Estimation has been performed:');
                    obj.Params.Coreg.out.regfile = reg;
                    
                    new = [outpath1 obj.Params.Coreg.in.prefix b1 c1];
                    obj.Params.Coreg.out.regimage = new;
                    
                    if exist(new,'file')==0
                        wasRun = true;
                        if obj.Params.Coreg.in.reslice
                            cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' source ' --reg ' reg ' --fstarg --o ' new];
                        else
                            cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' source ' --reg ' reg ' --fstarg --no-resample --o ' new];
                        end
                        runFS(cmd,obj.fsdir);
                    end

                    disp('coregistration with bbregister is complete');
                case 'iuw'    
                    disp('Perform this step with the IUW module');
                    return;
                otherwise
                    disp('unknown coregistration style option');
                    return
            end
                        
            if obj.Params.Coreg.in.getLabels == 1;
                obj.proc_get_fs_labels;
            end
            
            
            obj.UpdateHist(obj.Params.Coreg,'proc_coreg',new ,wasRun);
            fprintf('\n');
        end
        
        function obj = proc_apply_coreg(obj,fn,reslice)
            fprintf('\n%s\n', 'PERFORMING COREGISTRATION');
            wasRun = false;
                                    
            if ischar(fn)
                fn = cellstr(fn);
            end
            
            if nargin < 3
                reslice = obj.Params.Coreg.in.reslice;
            else
                obj.Params.Coreg.in.reslice = reslice;
            end
            
            [a b c] = fileparts(fn{1});
            if isempty(a); a = pwd; end;
            outpath1 = obj.getPath(a,obj.Params.Coreg.in.movefiles);
            outpath2 = obj.getPath(a,obj.Params.Coreg.in.movecrfiles);
            
            switch lower(obj.Params.Coreg.in.style)
                case 'spm'
                    disp('    USING SPM:');
                    load([outpath2 'CoReg.mat']);
                    
                    for ii = 1:numel(fn)
                        [a b c] = fileparts(fn{ii});
                        if isempty(a); a = pwd; end
                        nfn{ii,1} = [outpath1 obj.Params.Coreg.in.prefix b c];
                        if exist(nfn{ii},'file')==0
                            wasRun=true;
                            copyfile(fn{ii},nfn{ii});
                            MM = spm_get_space(nfn{ii});
                            spm_get_space(nfn{ii}, M\MM);
                            if exist([outpath1 obj.Params.Coreg.in.prefix b '.mat'])
                                delete([outpath1 obj.Params.Coreg.in.prefix b '.mat'])
                            end
                            
                            if obj.Params.Coreg.in.reslice
                                h1 = spm_vol(obj.Params.Coreg.in.target);
                                h2 = spm_vol(nfn{ii});
                                m = resizeVol2(h2,h1,obj.Params.Coreg.in.resampopt);
                                h1.fname = h2.fname;
                                h1.dt = h2.dt;
                                spm_write_vol(h1,m);
                            end
                        end
                        disp('Coregistration has been applied to the run:');
                    end
                    
                case 'bbreg'
                    disp('    USING BBREGISTER:');
                    reg = obj.Params.Coreg.out.regfile;
                    if isempty(reg)
                        error('no coreg file was specified');
                    end
                    
                    for ii = 1:numel(fn)
                        [a1 b1 c1] = fileparts(fn{ii});
                        if isempty(a1); a1 = pwd; end;
                        nfn{ii,1} = [outpath1 obj.Params.Coreg.in.prefix b1 c1];
                        if exist(nfn{ii,1},'file')==0
                            wasRun = true;
                            if obj.Params.Coreg.in.reslice
                                cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' fn{ii} ' --reg ' reg ' --fstarg --o ' nfn{ii}];
                            else
                                cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' fn{ii} ' --reg ' reg ' --fstarg --no-resample --o ' nfn{ii}];
                            end
                            runFS(cmd,obj.fsdir);
                        end
                        disp('Coregistration has been applied to the run:');
                    end

                    disp('coregistration with bbregister is complete');
                case 'iuw'    
                    disp('Perform this step with the IUW module');
                    return;
                otherwise
                    disp('unknown coregistration style option');
                    return
            end
            
            obj.Params.Coreg.out.otherin = fn;
            obj.Params.Coreg.out.otherout = nfn;

            obj.UpdateHist(obj.Params.Coreg,'proc_apply_coreg',nfn{1} ,wasRun);
            fprintf('\n');
        end
         
        function obj = proc_get_fs_labels(obj)
            fprintf('\n%s\n', 'FETCHING FS LABELS:');
            if strcmpi(obj.Params.FS_Params.regopt,'header');
                reg = ' --regheader ';
            else
                reg = [' --reg ' obj.Params.Coreg.out.regfile ' '];
            end
            
            if strcmpi(obj.Params.Coreg.in.style,'iuw')
                T1 = obj.Params.spmT1_Proc.in.t1;
                source = obj.Params.Implicit_Unwarp.out.newmean;
                [a b c] = fileparts(source);
                outpath = obj.getPath(a,obj.Params.Coreg.in.movecrfiles); 
                
               
                out1 = [outpath obj.Params.FS_Params.surf_out '_lh.nii'];
                cmd1 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' source ' ' reg T1 ' --hemi lh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out1 ];
                
                out2 = [outpath obj.Params.FS_Params.surf_out '_rh.nii'];
                cmd2 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' source ' ' reg T1 ' --hemi rh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out2 ];

                out3 = [outpath obj.Params.FS_Params.vol_out '.nii'];
                cmd3 = ['mri_label2vol --seg  '  obj.fsdir filesep obj.fsubj '/mri/' obj.Params.FS_Params.vol_labels ' --temp ' source ' ' reg T1 ' --fillthresh ' obj.Params.FS_Params.fillthresh ' --o '  out3];
            end
                        
            if strcmpi(obj.Params.Coreg.in.style,'bbreg')
                T1 = [obj.fsdir filesep obj.fsubj filesep 'mri' filesep 'T1.mgz'];
                [a b c] = fileparts(obj.Params.Coreg.in.source);
                outpath = obj.getPath(a,obj.Params.Coreg.in.movecrfiles);
                
                out1 = [outpath obj.Params.FS_Params.surf_out '_lh.nii'];
                cmd1 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --hemi lh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out1 ];
                
                out2 = [outpath obj.Params.FS_Params.surf_out '_rh.nii'];
                cmd2 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --hemi rh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out2 ];

                out3 = [outpath obj.Params.FS_Params.vol_out '.nii'];
                cmd3 = ['mri_label2vol --seg  '  obj.fsdir filesep obj.fsubj '/mri/' obj.Params.FS_Params.vol_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --fillthresh ' obj.Params.FS_Params.fillthresh ' --o '  out3];
            end            
            
            if strcmpi(obj.Params.Coreg.in.style,'spm')
                T1 = obj.Params.Coreg.in.target;
                [a b c] = fileparts(obj.Params.Coreg.in.source);
                outpath = obj.getPath(a,obj.Params.Coreg.in.movecrfiles);
                
                out1 = [outpath obj.Params.FS_Params.surf_out '_lh.nii'];
                cmd1 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --hemi lh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out1 ];
                
                out2 = [outpath obj.Params.FS_Params.surf_out '_rh.nii'];
                cmd2 = ['mri_label2vol --subject  '  obj.fsubj ' --annot ' obj.Params.FS_Params.surf_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --hemi rh --fillthresh ' obj.Params.FS_Params.fillthresh ' --proj frac ' obj.Params.FS_Params.proj_frac ' --o ' out2 ];

                out3 = [outpath obj.Params.FS_Params.vol_out '.nii'];
                cmd3 = ['mri_label2vol --seg  '  obj.fsdir filesep obj.fsubj '/mri/' obj.Params.FS_Params.vol_labels ' --temp ' obj.Params.Coreg.out.regimage reg T1 ' --fillthresh ' obj.Params.FS_Params.fillthresh ' --o '  out3];                
            end
            
            if exist(out1,'file')==0
                runFS(cmd1,obj.fsdir);
                Shell(out1,23);
            end
            
            if exist(out2,'file')==0
                runFS(cmd2,obj.fsdir);
                Shell(out2,23);
            end
            
            if exist(out3,'file')==0
                runFS(cmd3,obj.fsdir);
                Shell(out3,23);
            end
            
            obj.Params.Coreg.out.labels = {out1; out2; out3};            
            disp('FS labels have been mapped to functional space');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function obj = proc_fs_compute_vfm(obj)
            fprintf('\n%s\n', 'COMPUTE FS VOLUME FRACTION MAPS:');
            wasRun = false;
            
            if isempty(obj.fsubj)
                disp('No FS subject has been specified.');
                return
            end
            
            if isempty(obj.P_Coreg.source);
                disp('No source image has been specified');
                return
            end
            
            if isempty(obj.P_Coreg.regmean);
                disp('No registered image has been specified');
                return
            end
            
            [a b c] = fileparts(obj.P_Coreg.source);
            outpath = obj.getPath(a,obj.P_Coreg.movefiles);
            
            if exist([outpath 'areg.dat'],'file')==0
                cmd = ['tkregister2 --mov ' obj.P_Coreg.regmean ' --s ' obj.fsubj ' --regheader --reg ' outpath 'areg.dat --ltaout ' a filesep 'areg.lta --noedit'];
                runFS(cmd,obj.fsdir);
            end
            
            if exist([outpath 'vf.gm.nii'],'file')==0
                wasRun = true;
                cmd = ['mri_compute_volume_fractions ' outpath 'areg.dat ' obj.P_Coreg.regmean ' ' outpath 'vf;'];
                runFS(cmd,obj.fsdir);
                
                runFS(['mri_convert ' outpath 'vf.gm.mgz ' outpath 'vf.gm.nii'],obj.fsdir)
                runFS(['mri_convert ' outpath 'vf.wm.mgz ' outpath 'vf.wm.nii'],obj.fsdir)
                runFS(['mri_convert ' outpath 'vf.csf.mgz ' outpath 'vf.csf.nii'],obj.fsdir)
                
                m1 = (openIMG([outpath 'vf.gm.nii'])+openIMG([outpath 'vf.wm.nii']))/2;
                i1 = find(m1==0);
                [m2 h] = openIMG([outpath 'vf.csf.nii']);
                m2(i1)=0;
                spm_write_vol(h,m2);
                
                delete([outpath 'vf*.mgz']);
                % delete([outpath 'areg.dat']);
            end
            
            obj.P_Coreg.bbreg.vfmaps = { ...
                [outpath 'vf.gm.nii'];
                [outpath 'vf.wm.nii'];
                [outpath 'vf.csf.nii'];
            };
            disp('volume fraction maps have been computed');
            
            obj.UpdateHist(obj.P_FS_Params,'proc_fs_compute_vfm', [outpath 'vf.gm.nii'],wasRun);
        end
        
        function obj = proc_filter(obj,fn)
            fprintf('\n\n%s\n', 'FILTERING DATA:');
            %%% MODULE IS COMPLETE
            wasRun = false;
            
            if nargin<2
                fn = obj.Params.Filter.in.fn;
            else
                obj.Params.Filter.in.fn = fn;
            end
            
            nfn = [];
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.Params.Filter.in.movefiles);                
                nfn{jj,1} = [outpath obj.Params.Filter.in.prefix b c];
            end
            
            for ii = 1:length(fn)
                if exist(nfn{ii},'file')==0
                    wasRun = true;
                    V = spm_vol([fn{ii} ',1']);
                    M = FastRead(fn{ii});
                    
                    if obj.Params.Filter.in.ampscale
                        disp('scaling by mean bold');
                        base = nanmean(M(:)); 
                        sc = base+abs(mean(M));
                        M = (M./repmat(sc,size(M,1),1))*100;
                    end
                    M(~isfinite(M))=0;
                    
                    mn = nanmean(M);
                    
                    if obj.Params.Filter.in.detrend == 1
                        disp('applying linear detrend');
                        M = detrend(M);
                    end
                    
                    if isnan(obj.Params.Filter.in.lowcut)
                        disp('applying high-pass filter');
                        M = ft_preproc_highpassfilter(M',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                    elseif isnan(obj.Params.Filter.in.highcut)
                        disp('applying low-pass filter');
                        M = ft_preproc_lowpassfilter(M',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                    else
                        disp('applying band-pass filter');
                        M = ft_preproc_bandpassfilter(M',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                    end
                    
                    
                    
                    M = reshape(M', V.dim(1), V.dim(2), V.dim(3), size(M,1));
                    V.fname = nfn{ii};
                    if V.dt(1)<16
                        V.dt(1) = 16;
                    end
                    
                    for jj = 1:size(M,4)
                        V.n = [jj 1];
                        spm_write_vol(V,M(:,:,:,jj));
                    end
                end
                disp('temporal filtering of run is complete.');
            end
            
            obj.Params.Filter.in.fn = fn;
            obj.Params.Filter.out.fn = nfn;
            obj.lastFN = nfn;

            obj.UpdateHist(obj.Params.Filter,'proc_filter',nfn{end},wasRun);
            fprintf('\n');
        end

        function obj = proc_compute_physio_regs(obj,fn)
            fprintf('\n%s\n', 'CREATE PHYSIO NUISANCE REGRESSOR FILES:');
            wasRun=false;
            
            if nargin<2
                fn = obj.Params.ComputePhysioRegs.in.fn;
            else
                obj.Params.ComputePhysioRegs.in.fn = fn;
            end
            
            if ischar(fn)
                fn = cellstr(fn);
            end
            
            switch lower(obj.Params.ComputePhysioRegs.in.type)
                case 'masks'
                    disp('Use predifined masks:');
                    obj.Params.ComputePhysioRegs.out.regs = [];
                    obj.Params.ComputePhysioRegs.out.resfiles = [];
                    for ii = 1:numel(fn)
                        [a b c] = fileparts(fn{ii});
                        outpath = obj.getPath(a,obj.Params.ComputePhysioRegs.in.movefiles);
                        obj.Params.ComputePhysioRegs.out.resfiles{ii} = [outpath b '_' obj.Params.ComputePhysioRegs.in.filename];
                        if exist(obj.Params.ComputePhysioRegs.out.resfiles{ii},'file')>0
                            dat = load(obj.Params.ComputePhysioRegs.out.resfiles{ii});
                            obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                            continue
                        end
                        
                        wasRun = true;
                        
                        V = spm_vol([fn{ii} ',1']);
                        m = FastRead(fn{ii});
                        
                        dat = [];
                        for jj = 1:numel(obj.Params.ComputePhysioRegs.in.masks)
                            h1 = spm_vol(obj.Params.ComputePhysioRegs.in.masks{jj});
                            
                            map = resizeVol2(h1,V,obj.Params.ComputePhysioRegs.in.resample);
                            map = map(:)';
                            
                            if obj.Params.ComputePhysioRegs.in.weighted == 1
                              mm = m.*repmat(map,size(m,1),1);
                              sc = sum(map);
                            else
                                %%% might need to allow for different thresholds for different maps
                                if numel(obj.Params.ComputePhysioRegs.in.threshold)==1
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold);
                                else
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold(jj));
                                end
                                mm = m.*repmat(ind(:)',size(m,1),1);
                                sc = sum(ind);
                            end
                                                        
                            if obj.Params.ComputePhysioRegs.in.prinComps
                                if numel(obj.Params.ComputePhysioRegs.in.nPC)==1
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC,1000);
                                else
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC(jj),1000);
                                end
                                dat = [dat PCs];
                            else
                                dat(:,end+1) = sum(mm,2)./sc;
                            end
                        end
                        
                        obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                        save(obj.Params.ComputePhysioRegs.out.resfiles{ii},'dat','-ascii');
                    end
                case 'fsvf'
                    disp('Use FreeSurfer Volume Fraction maps:');
                    obj.Params.ComputePhysioRegs.out.regs = [];
                    obj.Params.ComputePhysioRegs.out.resfiles = [];
                    
                    for ii = 1:numel(fn)         
                        [a b c] = fileparts(fn{ii});
                        outpath = obj.getPath(a,obj.Params.ComputePhysioRegs.in.movefiles);
                        obj.Params.ComputePhysioRegs.out.resfiles{ii} = [outpath b '_' obj.Params.ComputePhysioRegs.in.filename];
                        if exist(obj.Params.ComputePhysioRegs.out.resfiles{ii},'file')>0
                            dat = load(obj.Params.ComputePhysioRegs.out.resfiles{ii});
                            obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                            continue
                        end
            
                        wasRun = true;
                        
                        V = spm_vol([fn{ii} ',1']);
                        m = FastRead(fn{ii});
                        
                        dat = [];
                        for jj = obj.Params.ComputePhysioRegs.in.whichparts(:)'                            
                            h1 = spm_vol(obj.P_Coreg.bbreg.vfmaps{jj});
                            if ~all(all(round(h1.mat*10000)==round(V.mat*10000)))
                                error('Estimate TPM maps have a different header than the input files');
                            end
                            
                            map = resizeVol2(h1,V,obj.Params.ComputePhysioRegs.in.resample);
                            map = map(:)';
                            
                            if obj.Params.ComputePhysioRegs.in.weighted == 1
                              mm = m.*repmat(map,size(m,1),1);
                              sc = sum(map);
                            else
                                %%% might need to allow for difference thresholds for different maps
                                if numel(obj.Params.ComputePhysioRegs.in.threshold)==1
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold);
                                else
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold(jj));
                                end
                                mm = m.*repmat(ind(:)',size(m,1),1);
                                sc = sum(ind);
                            end
                                                        
                            if obj.Params.ComputePhysioRegs.in.prinComps
                                if numel(obj.Params.ComputePhysioRegs.in.nPC)==1
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC,1000);
                                else
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC(jj),1000);
                                end
                                dat = [dat PCs];
                            else
                                dat(:,end+1) = sum(mm,2)./sc;
                            end
                        end

                        obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                        save(obj.Params.ComputePhysioRegs.out.resfiles{ii},'dat','-ascii');
                    end
                case 'directtpm' 
                    disp('using the DirectTPM approach based on the direct normalization');
                    obj.Params.ComputePhysioRegs.out.regs = [];
                    obj.Params.ComputePhysioRegs.out.resfiles = [];
                    
                    for ii = 1:numel(fn)         
                        [a b c] = fileparts(fn{ii});
                        outpath = obj.getPath(a,obj.Params.ComputePhysioRegs.in.movefiles);
                        obj.Params.ComputePhysioRegs.out.resfiles{ii} = [outpath b '_' obj.Params.ComputePhysioRegs.in.filename];
                        if exist(obj.Params.ComputePhysioRegs.out.resfiles{ii},'file')>0
                            dat = load(obj.Params.ComputePhysioRegs.out.resfiles{ii});
                            obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                            continue
                        end
            
                        wasRun = true;
                        
                        V = spm_vol([fn{ii} ',1']);
                        m = FastRead(fn{ii});
                        
                        dat = [];
                        for jj = obj.Params.ComputePhysioRegs.in.whichparts(:)'
                            h1 = spm_vol(obj.Params.NormNew.out.estTPM{jj});
                            if ~all(all(round(h1.mat*10000)==round(V.mat*10000)))
                                error('Estimate TPM maps have a different header than the input files');
                            end
                            
                            map = FastRead(obj.Params.NormNew.out.estTPM{jj});

                            if obj.Params.ComputePhysioRegs.in.weighted == 1
                              mm = m.*repmat(map,size(m,1),1);
                              sc = sum(map);
                            else
                                %%% might need to allow for difference thresholds for different maps
                                if numel(obj.Params.ComputePhysioRegs.in.threshold)==1
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold);
                                else
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold(jj));
                                end
                                mm = m.*repmat(ind(:)',size(m,1),1);
                                sc = sum(ind);
                            end
                                                        
                            if obj.Params.ComputePhysioRegs.in.prinComps
                                %T = ft_preproc_bandpassfilter(T',1/3,[0.01 0.10],4,'but','twopass' ,'no')';
                                if numel(obj.Params.ComputePhysioRegs.in.nPC)==1
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC,1000);
                                else
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC(jj),1000);
                                end
                                dat = [dat PCs];
                            else
                                dat(:,end+1) = sum(mm,2)./sc;
                            end
                        end

                        obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                        save(obj.Params.ComputePhysioRegs.out.resfiles{ii},'dat','-ascii');
                    end   
                case 'indirecttpm'
                    disp('using the IndirectTPM approach based on the SPM T1 segmentation');
                    obj.Params.ComputePhysioRegs.out.regs = [];
                    obj.Params.ComputePhysioRegs.out.resfiles = [];
                    
                    for ii = 1:numel(fn)
                        [a b c] = fileparts(fn{ii});
                        outpath = obj.getPath(a,obj.Params.ComputePhysioRegs.in.movefiles);
                        obj.Params.ComputePhysioRegs.out.resfiles{ii} = [outpath b '_' obj.Params.ComputePhysioRegs.in.filename];
                        if exist(obj.Params.ComputePhysioRegs.out.resfiles{ii},'file')>0
                            dat = load(obj.Params.ComputePhysioRegs.out.resfiles{ii});
                            obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                            continue
                        end
                        
                        wasRun = true;
                        
                        V = spm_vol([fn{ii} ',1']);
                        m = FastRead(fn{ii});
                        
                        dat = [];
                        for jj = obj.Params.ComputePhysioRegs.in.whichparts(:)'
                            
                            h1 = spm_vol(obj.Params.spmT1_Proc.out.estTPM{jj});
                            
                            map = resizeVol2(h1,V,obj.Params.ComputePhysioRegs.in.resample);
                            map = map(:)';
                            
                            if obj.Params.ComputePhysioRegs.in.weighted == 1
                                mm = m.*repmat(map,size(m,1),1);
                                sc = sum(map);
                            else
                                %%% might need to allow for difference thresholds for different maps
                                if numel(obj.Params.ComputePhysioRegs.in.threshold)==1
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold);
                                else
                                    ind = (map>obj.Params.ComputePhysioRegs.in.threshold(jj));
                                end
                                mm = m.*repmat(ind(:)',size(m,1),1);
                                sc = sum(ind);
                            end
                            
                            if obj.Params.ComputePhysioRegs.in.prinComps
                                if numel(obj.Params.ComputePhysioRegs.in.nPC)==1
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC,1000);
                                else
                                    PCs = IterPCA(demean(mm),obj.Params.ComputePhysioRegs.in.nPC(jj),1000);
                                end
                                dat = [dat PCs];
                            else
                                dat(:,end+1) = sum(mm,2)./sc;
                            end
                        end
                        obj.Params.ComputePhysioRegs.out.regs{ii} = dat;
                        save(obj.Params.ComputePhysioRegs.out.resfiles{ii},'dat','-ascii');
                    end
                otherwise
                    error('type option was not recognized');
            end
            obj.Params.ComputePhysioRegs.in.fn = fn;
            
            obj.UpdateHist(obj.Params.ComputePhysioRegs,'proc_compute_physio_regs',obj.Params.ComputePhysioRegs.out.resfiles{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_map_to_fsaverage(obj,fn)
            fprintf('\n%s\n', 'MAP TO FSAVERAGE:');
            wasRun = false;
                        
            if ischar(fn)
                fn = cellstr(fn);
            end
            
            
            if strcmpi(obj.Params.FS_Params.regopt,'header');
                reg = ' --regheader ';
            else
                reg = [' --reg ' obj.Params.Map_to_Surface.in.regfile ' '];
            end
            
            nfn = [];
            for ii = 1:numel(fn)
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.Params.Map_To_Surface.in.movefiles);

                nfn{ii,1} = [outpath 'ss' obj.Params.FS_Params.surf_ss '_' b '_lh' c];
                nfn{ii,2} = [outpath 'ss' obj.Params.FS_Params.surf_ss '_' b '_rh' c];
                
                if exist(nfn{ii,1}, 'file')>0
                    continue;
                end
                wasRun = true;
                
                C1 = ['mri_vol2surf ' ...
                    '--mov ' fn{ii} ' ' ...
                    reg obj.fsubj ' '  ...
                    '--trgsubject ' obj.Params.FS_Params.fsaverage ' ' ...
                    '--projfrac '   obj.Params.FS_Params.projfrac ' ' ...
                    '--surf-fwhm '  obj.Params.FS_Params.surf_ss ' ' ...
                    '--hemi lh ' ...
                    '--out ' nfn{ii,1} ';'];
                
                C2 = ['mri_vol2surf ' ...
                    '--mov ' fn{ii} ' ' ...
                    reg obj.fsubj ' '  ...
                    '--trgsubject ' obj.Params.FS_Params.fsaverage ' ' ...
                    '--projfrac '   obj.Params.FS_Params.projfrac ' ' ...
                    '--surf-fwhm '  obj.Params.FS_Params.surf_ss ' ' ...
                    '--hemi rh ' ...
                    '--out ' nfn{ii,2} ';'];
                
                runFS([C1 C2],obj.fsdir);
            end
            
            obj.Params.Map_To_Surface.in.fn = fn;
            obj.Params.Map_To_Surface.out.fn = nfn;
            
            obj.UpdateHist( obj.Params.Map_To_Surface,'proc_map_to_fsaverage',nfn{end},wasRun);
            fprintf('\n');
        end
       
        function obj = proc_implict_unwarping(obj,source,fn)
            fprintf('\n%s\n', 'WARP DATA TO ITS OWN T1 IMAGE:');
            wasRun = false;
            
            bb = world_bb(spm_vol(obj.Params.spmT1_Proc.in.t1));
            
            if nargin<2
                source = obj.Params.Implicit_Unwarp.in.source;
            else
                obj.Params.Implicit_Unwarp.in.source = source;
            end
            
            
            if nargin<3
                fn = [];
                fn = obj.Params.Implicit_Unwarp.in.fn;
            else
                if ischar(fn)
                    fn = cellstr(fn);
                end
                
                obj.Params.Implicit_Unwarp.in.fn = fn;
            end
            
            [a b c] = fileparts(source);
            outpath = obj.getPath(a,obj.Params.Implicit_Unwarp.in.movefiles);
           
            
            localMean = [outpath 'mean.nii'];
            %%% Copy over the mean EPI
            if exist([outpath 'mean.nii'],'file')==0
                [tm th] = openIMG(source);
                th.fname = localMean;
                tm(isnan(tm))=0; spm_write_vol(th,tm);
            end
            disp('copy mean image');

            %%% Create the localTPM.nii maps
            if ~exist([outpath '/localTPM.nii']);
                ll = obj.Params.spmT1_Proc.out.estTPM;
                for ii = 1:numel(ll)
                    spm_smooth(ll{ii},[outpath 'localTPM.nii,' num2str(ii)],[repmat(obj.Params.Implicit_Unwarp.in.tpmSmooth,1,3)],16);
                end
            end
            disp('create local tpm');
            
            %%% Perform the Affine registration
            if exist([outpath 'Affine.mat'],'file')==0
                % tpm  = spm_load_priors8(char(obj.P_spmT1_Proc.estTPM));
                tpm  = spm_load_priors8([outpath 'localTPM.nii']);
                
                in.P = localMean;
                in.samp = 3;
                in.fwhm1 = 3;
                in.fwhm2 = 0;
                in.tpm = tpm;
                in.M = [];
                in.regtype='rigid';
                
                [Affine,h] = spm_maff8(in.P,in.samp,in.fwhm1,in.tpm,in.M,in.regtype);
                in.M = Affine;
                [Affine,h] = spm_maff8(in.P,in.samp,in.fwhm2,in.tpm,in.M,in.regtype);
                save([outpath 'Affine.mat'],'Affine');
            else
                load([outpath 'Affine.mat'],'Affine');
            end
            disp('compute affine transformation');
            
            %%% Compute the nonlinear registration
            if exist([outpath 'Norm_Seg8.mat'],'file')==0
                tpm  = spm_load_priors8([outpath 'localTPM.nii']);
                
                P.nn.NormPars.image = spm_vol(localMean);
                P.nn.NormPars.fwhm = 0;
                P.nn.NormPars.biasreg = 0.0001;
                P.nn.NormPars.biasfwhm = 60;
                P.nn.NormPars.tpm = tpm;
                P.nn.NormPars.lkp = [];
                P.nn.NormPars.reg = [0 0.001 0.5 0.05 0.2];
                %P.nn.NormPars.reg = [1 1 1 1 1];
                %P.nn.NormPars.reg = [0.001 0.001 0.001 0.001 0.001];
                P.nn.NormPars.samp = 3;
                P.nn.NormPars.Affine = Affine;
                
                results = spm_preproc8(P.nn.NormPars);
                
                save([outpath 'Norm_Seg8.mat'], 'results');
            else
                load([outpath 'Norm_Seg8.mat'], 'results');
            end
            disp('compute normalization');
            
            %%%
            if exist([outpath  'y_mean.nii'],'file')==0
                
                w = [0 0 0 0;
                     0 0 0 0;
                     0 0 0 0;
                     0 0 0 0;
                     0 0 0 0;
                     0 0 0 0];
                [cls,M1] = spm_preproc_write8(results,w,[1 1],[1 1],0,1,bb,obj.vox(1));
            end
            obj.Params.Implicit_Unwarp.out.regfile =  [outpath  'y_mean.nii'];
            obj.Params.Implicit_Unwarp.out.iregfile = [outpath 'iy_mean.nii'];
            
            
            %%% set defs
            defs.comp{1}.def = {obj.Params.Implicit_Unwarp.out.regfile};
            defs.comp{2}.idbbvox.vox = obj.vox;
            defs.comp{2}.idbbvox.bb = bb;
            
            defs.out{1}.pull.fnames = [];
            defs.out{1}.pull.savedir.savesrc = 1;
            defs.out{1}.pull.interp=obj.interporder;
            defs.out{1}.pull.mask=1;
            defs.out{1}.pull.fwhm=[0 0 0];
            
            if exist([outpath obj.Params.Implicit_Unwarp.in.prefix 'mean.nii'],'file')==0
                defs.out{1}.pull.fnames = {[outpath 'mean.nii']};
               
                spm_deformations(defs);
                movefile([outpath 'wmean.nii'],[outpath obj.Params.Implicit_Unwarp.in.prefix 'mean.nii']);
            end
            obj.Params.Implicit_Unwarp.out.defs = defs;
            %%%
            nfn = [];
            for ii = 1:numel(fn)
                [a b c] = fileparts(fn{ii});
                if exist([outpath obj.Params.Implicit_Unwarp.in.prefix b c],'file')==0
                    defs.out{1}.pull.fnames = fn(ii);
                    spm_deformations(defs);
                    
                    movefile([a filesep 'w' b c],[outpath obj.Params.Implicit_Unwarp.in.prefix b c]);
                    nfn{ii,1} = [outpath obj.Params.Implicit_Unwarp.in.prefix b c];
                else
                    nfn{ii,1} = [outpath obj.Params.Implicit_Unwarp.in.prefix b c];
                end
            end
            
            
            %%%
            disp('apply normalization');
            obj.Params.Implicit_Unwarp.out.fn = nfn;
            obj.Params.Implicit_Unwarp.out.newmean = [outpath obj.Params.Implicit_Unwarp.in.prefix 'mean.nii'];

            obj.UpdateHist(obj.Params.Implicit_Unwarp,'proc_implict_unwarping',obj.Params.Implicit_Unwarp.out.regfile,wasRun);
        end
        
        function obj = proc_regress_clean(obj,fn)
            fprintf('\n\n%s\n', 'CLEANING DATA:');
            wasRun = false;
            
            if ischar(fn)
                fn = cellstr(fn);
            end

            nfn = [];
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.Params.CleanData.in.movefiles);
                if exist(outpath,'dir')==0; mkdir(outpath); end
                
                nfn{jj,1} = [outpath obj.Params.CleanData.in.prefix b c];
            end
            
            obj.Params.CleanData.out.REGS = [];
            REGS = [];
            for ii = 1:length(fn)
                if exist(nfn{ii},'file')==0
                    wasRun = true;
                    V = spm_vol([fn{ii} ',1']);
                    M = FastRead(fn{ii});
                                        
                    if obj.Params.CleanData.in.filter
                        if obj.Params.Filter.in.detrend == 1
                            disp('applying linear detrend to the data');
                            M = detrend(M);
                        end
                        
                        if isnan(obj.Params.Filter.in.lowcut)
                            disp('applying high-pass filter');
                            M = ft_preproc_highpassfilter(M',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        elseif isnan(obj.Params.Filter.in.highcut)
                            disp('applying low-pass filter');
                            M = ft_preproc_lowpassfilter(M',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        else
                            disp('applying band-pass filter');
                            M = ft_preproc_bandpassfilter(M',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        end
                        
%                         if obj.Params.Filter.in.highcut == 0
%                             disp('applying high-pass filter to the data');
%                             M = ft_preproc_highpassfilter(M',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         elseif obj.Params.Filter.in.lowcut == 0
%                             disp('applying low-pass filter to the data');
%                             M = ft_preproc_lowpassfilter(M',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         else
%                             disp('applying band-pass filter to the data');
%                             M = ft_preproc_bandpassfilter(M',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         end
                    end
                    
                    regs = [];
                    if obj.Params.CleanData.in.motion
                        regs = [regs load(obj.Params.Realign.out.realigpars{ii})];
                    end
                    
                    if obj.Params.CleanData.in.physio
                        regs = [regs obj.Params.ComputePhysioRegs.out.regs{ii}];
                    end
                    
                    if obj.Params.CleanData.in.other    
                        regs = [regs load(obj.Params.CleanData.in.otherFiles{ii})];
                    end
                    
                    if obj.Params.CleanData.in.deriv
                        regs = [regs [zeros(1,size(regs,2)); diff(regs)]];
                    end
                    
                    if obj.Params.CleanData.in.square
                        regs = [regs regs.^2];
                    end
                    
                    
                    if obj.Params.CleanData.in.filter
                        if obj.Params.Filter.in.detrend == 1
                            disp('applying linear detrend to the regressors');
                            regs = detrend(regs);
                        end
                        
                        if isnan(obj.Params.Filter.in.lowcut)
                            disp('applying high-pass filter to the regressors');
                            regs = ft_preproc_highpassfilter(regs',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        elseif isnan(obj.Params.Filter.in.highcut)
                            disp('applying low-pass filter to the regressors');
                            regs = ft_preproc_lowpassfilter(regs',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        else
                            disp('applying band-pass filter to the regressors');
                            regs = ft_preproc_bandpassfilter(regs',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        end
                        
%                         if obj.Params.Filter.in.highcut == 0
%                             disp('applying high-pass filter to the regressors');
%                             regs = ft_preproc_highpassfilter(regs',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         elseif obj.Params.Filter.in.lowcut == 0
%                             disp('applying low-pass filter to the regressors');
%                             regs = ft_preproc_lowpassfilter(regs',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         else
%                             disp('applying band-pass filter to the regressors');
%                             regs = ft_preproc_bandpassfilter(regs',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
%                         end
                    end
                    
                    disp(['there are ' num2str(size(regs,2)) ' regressors in the set']);
                    
                    if obj.Params.CleanData.in.reduce
                        disp('performing data reduction on the regressor set ...');
                        regs = zscore(regs);
                        
                        cm = (regs'*regs)./(size(regs,1)-1);
                        [U,sigma,RM] = svd(cm);
                        
                        
                        tmp = abs([max(RM); min(RM)]);
                        [trash,i1] = max(tmp); i1(i1==2)=-1;
                        for kk = 1:size(RM,2); RM(:,kk) = RM(:,kk)*i1(kk); end
                        
                        
                        E = diag(sigma);
                        PCs = regs*RM;
                        
                        %%% must account for at least 3% of the variance
                        %thresh = (size(regs,2)/100)*3;
                        %regs = PCs(:,find(E>thresh));
                        
                        kk = find(cumsum(E)<(.9*numel(E))); %kk(end)
                        regs = PCs(:,kk);
                        
                        %[freq, amp, coeff] = power_spec(PCs(:,1), 1/3);
                        %[freq, amp, coeff] = power_spec(PCs(:,1), 1/3);
                    end
                    disp(['there are ' num2str(size(regs,2)) ' regressors in the set']);
                    obj.Params.CleanData.out.REGS{ii} = regs;
                    [a b c] = fileparts(fn{ii});
                    save([outpath b '_regs.txt'],'regs','-ascii');
                    
                    %%%
                    pred = regs*(pinv(regs)*M);
                    vol = zeros(V.dim);
                    vol(:) = SumOfSquares(pred)./SumOfSquares(M);
                    V.fname = [outpath b '_R2reg.nii'];
                    V.dt = [16 0];
                    spm_write_vol(V,vol);
                    %%%
                    
                    if size(regs,2)>0
                        M = crtlFor(M,regs);
                    end
                    
                    M = reshape(M', V.dim(1), V.dim(2), V.dim(3), size(M,1));
                    V.fname = nfn{ii};
                    if V.dt(1)<16
                        V.dt(1) = 16;
                    end
                    
                    for jj = 1:size(M,4)
                        V.n = [jj 1];
                        spm_write_vol(V,M(:,:,:,jj));
                    end
                end
                disp('temporal filtering of run is complete.');
            end
            
            if ~exist([[outpath filesep 'R2map.nii']],'file')
                obj.proc_average(dir_wfp([outpath filesep '*R2reg.nii']),'R2map.nii');
            end
            
            obj.Params.CleanData.in.fn = fn;
            obj.Params.CleanData.out.fn = nfn;
            obj.lastFN = nfn;

            obj.UpdateHist(obj.Params.CleanData,'proc_regress_clean',nfn{end},wasRun);
            fprintf('\n');
        end

        function obj = proc_regress_clean2(obj,fn)
            fprintf('\n\n%s\n', 'CLEANING DATA v2:');
            wasRun = false;
            %%% This section is to experiment with a unitary approach to
            %%% nuisance regression across time, that is nuisance variables
            %%% are regressed out once with a common set of parameter
            %%% estimates rather than independent parameter estimation for
            %%% each run.
            nfn = []; check = 0;
            for jj = 1:numel(fn)
                [a b c] = fileparts(fn{jj});
                outpath = obj.getPath(a,obj.Params.CleanData.in.movefiles);                
                nfn{jj,1} = [outpath obj.Params.CleanData.in.prefix b c];
                check = check + (exist(nfn{jj},'file')>0);
            end
            
            if check == numel(nfn)
                disp('processing is complete');
                obj.Params.CleanData.in.fn = fn;
                obj.Params.CleanData.out.fn = nfn;
                obj.lastFN = nfn;
                
                obj.UpdateHist(obj.Params.CleanData,'proc_regress_clean2',nfn{end},wasRun);
                return
            end
            
            %%%
            M = [];
            REGS = []; C = []; II = [];
            for ii = 1:numel(fn);
                m = FastRead(fn{ii});
                si1 = find(~isnan(sum(m)) & std(m)>0);
                
                regs = [];
                if obj.Params.CleanData.in.motion
                    regs = [regs load(obj.Params.Realign.out.realigpars{ii})];
                end
                    
                if obj.Params.CleanData.in.physio
                    regs = [regs obj.Params.ComputePhysioRegs.out.regs{ii}];
                end
                    
                if obj.Params.CleanData.in.other
                    regs = [regs load(obj.Params.CleanData.in.otherFiles{ii})];
                end
                
                if obj.Params.CleanData.in.deriv
                    regs = [regs [zeros(1,size(regs,2)); diff(regs)]];
                end
                
                if obj.Params.CleanData.in.square
                    regs = [regs regs.^2];
                end
                    
                if obj.Params.CleanData.in.filter
                    if obj.Params.Filter.in.detrend == 1
                        disp('applying linear detrend');
                        regs = detrend(regs);
                        m(:,si1) = detrend(m(:,si1));
                    end
                    if isnan(obj.Params.Filter.in.lowcut)
                        disp('applying high-pass filter');
                        regs = ft_preproc_highpassfilter(regs',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        %regs = ft_preproc_highpassfilter(regs',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        
                        m(:,si1) = ft_preproc_highpassfilter(m(:,si1)',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        %m(:,si1) = ft_preproc_highpassfilter(m(:,si1)',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        
                    elseif isnan(obj.Params.Filter.in.highcut)
                        disp('applying low-pass filter');
                        regs = ft_preproc_lowpassfilter(regs',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        %regs = ft_preproc_lowpassfilter(regs',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        
                        m(:,si1) = ft_preproc_lowpassfilter(m(:,si1)',1/obj.TR,obj.Params.Filter.in.lowcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        %m(:,si1) = ft_preproc_lowpassfilter(m(:,si1)',1/obj.TR,obj.Params.Filter.in.highcut,obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                    else
                        disp('applying band-pass filter');
                        regs = ft_preproc_bandpassfilter(regs',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        % regs = ft_preproc_bandpassfilter(regs',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        
                        m(:,si1) = ft_preproc_bandpassfilter(m(:,si1)',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                        %m(:,si1) = ft_preproc_bandpassfilter(m(:,si1)',1/obj.TR,[obj.Params.Filter.in.lowcut obj.Params.Filter.in.highcut],obj.Params.Filter.in.filterorder,'but','twopass' ,'no')';
                    end
                end
                
                II{ii} = ((size(REGS,1)+1): size(REGS,1)+size(regs,1))';
                REGS = [REGS; zscore(regs)];
                C = [C; zeros(size(regs,1),1)+ii];
                
                M = [M; demean(m)];
            end
            %%%
            disp(['there are ' num2str(size(REGS,2)) ' regressors in the set']);
            if obj.Params.CleanData.in.reduce
                regs = zscore(REGS);
                
                cm = (regs'*regs)./(size(regs,1)-1);
                [U,sigma,RM] = svd(cm);
                
                
                tmp = abs([max(RM); min(RM)]);
                [trash,i1] = max(tmp); i1(i1==2)=-1;
                for kk = 1:size(RM,2); RM(:,kk) = RM(:,kk)*i1(kk); end
                
                
                E = diag(sigma);
                PCs = regs*RM;
                
                % thresh = 1;
                % regs = PCs(:,find(E>thresh));
                
                kk = find(cumsum(E)<(.9*numel(E))); kk(end)
                regs = PCs(:,kk);
                
                
                %[freq, amp, coeff] = power_spec(PCs(:,1), 1/3);
                %[freq, amp, coeff] = power_spec(PCs(:,1), 1/3);
                
                %disp([sprintf('%3.0f',(sum(E>thresh)/numel(E))*100) '% of the variance was retained']);
                disp([sprintf('%3.0f',(sum(E(kk))/numel(E))*100) '% of the variance was retained']);
            else
                regs = REGS;
            end
            disp(['there are ' num2str(size(regs,2)) ' regressors in the set']);
            obj.Params.CleanData.out.REGS = regs;
            obj.Params.CleanData.out.indices = II;
            
            [a b c] = fileparts(fn{1});
            save([outpath b '_regs.txt'],'regs','-ascii');
            
            %%%
            h = spm_vol([fn{1} ',1']);
            vol = zeros(h.dim);
            
            pred = regs*(pinv(regs)*M);
            vol(:) = SumOfSquares(pred)./SumOfSquares(M);
            h.dt = [16 0];
            h.fname = [outpath 'R2map.nii'];
            spm_write_vol(h,vol);
            %%%
            M = crtlFor(M,regs);
            V = spm_vol([fn{1} ',1']);
            %%%
            for ii = 1:numel(II)
            
                mm = reshape(M(II{ii},:)', V.dim(1), V.dim(2), V.dim(3), numel(II{ii}));
                V.fname = nfn{ii};
                if V.dt(1)<16
                    V.dt(1) = 16;
                end
            
                for jj = 1:size(mm,4)
                    V.n = [jj 1];
                    spm_write_vol(V,mm(:,:,:,jj));
                end
            end
            
            obj.Params.CleanData.in.fn = fn;
            obj.Params.CleanData.out.fn = nfn;
            obj.lastFN = nfn;

            obj.UpdateHist(obj.Params.CleanData,'proc_regress_clean2',nfn{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_zscore(obj,fn)
            fprintf('\n\n%s\n', ['Z-Scoring Across the 4th dimension']);
            %%% MODULE IS COMPLETE
            wasRun = false;
            nfn = [];
            %%%
            for ii = 1:numel(fn);
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.P_ZScore.movefiles);
                if exist(outpath,'dir')==0
                    mkdir(outpath);
                end
                
                nfn{ii,1} = [outpath obj.P_ZScore.prefix b c];
                
                if exist(nfn{ii},'file')>0
                    disp('specified volumes have been dropped');
                    continue
                end
                
                wasRun = true;
                
                [M V] = openIMG(fn{ii});
                M = zscore(M,1,4);
                
                
                V = V(1);
                cc = 0;
                for jj = 1:size(M,4)
                    cc = cc+1;                    
                    V.fname = nfn{ii,1};
                    V.n = [cc 1];
                    spm_write_vol(V,M(:,:,:,jj));
                end
                disp('data has been statistically normalized');
            end
            obj.P_ZScore.in = fn;
            obj.P_ZScore.out = nfn;
            obj.lastFN = nfn;
            
            
            obj.UpdateHist(obj.P_ZScore,'proc_zscore',nfn{end},wasRun);
            fprintf('\n');
        end
        
        function obj = proc_average(obj,fn,new)
            if ischar(fn)
                fn = cellstr(fn);
            end
            
            [a b c] = fileparts(fn{1});
            if nargin == 2
                out  = [a filesep 'stdev_' b c];
            else
                out  = [a filesep new];
            end
            
            [M h] = openIMG(char(fn));
            h = h(1);
            h.dt = [16 0];
            h.fname = out;
            
            M = nanmean(M,4);
            spm_write_vol(h,M);
        end
        
        function obj = proc_stdev(obj,fn,new)
            if ischar(fn)
                fn = cellstr(fn);
            end
            
            [a b c] = fileparts(fn{1});
            if nargin == 2
                out  = [a filesep 'stdev_' b c];
            else
                out  = [a filesep new];
            end
            
            [M h] = openIMG(char(fn));
            h = h(1);
            h.dt = [16 0];
            h.fname = out;
            
            M = nanstd(M,1,4);
            spm_write_vol(h,M)
        end
        
        function obj = proc_wta_tpm(obj)
            fn = obj.Params.spmT1_Proc.out.estTPM;
            [a b c] = fileparts(fn{1});
            
            h = spm_vol(fn{1});
            h.dt = [4 0];
            h.pinfo(1) = 1;
            
            m = FastRead(fn);
            
            [tr i] = max(m);
            
            vol = zeros(h.dim);
            vol(:) = i;
            
            h.fname = [a filesep 'wta.nii'];
            spm_write_vol(h,vol);
        end
        
        function obj = proc_qa(obj)
            fprintf('\n%s\n', 'Examing Image Files and Creating Bad Volume Regressors:');
            
            if ~isempty(obj.Params.QA.out.reasons)
                disp('QA has already been preformed');
                return;
            end
            %%% gather information
            fn = obj.Params.SNR.in.fn;
            obj.Params.QA.out.reasons = table;
            for ii = 1:numel(fn)
                [a b c] = fileparts(fn{ii});
                outpath = obj.getPath(a,obj.Params.QA.in.movefiles);
                
                nnn = [outpath obj.Params.QA.in.logname '_' b '.txt'];
                
                snr = dlmread(obj.Params.SNR.out.report{ii},'\t',1,1);
                disp(snr(1:3));
                obj.Params.QA.out.SNR(ii,1) = snr(3);
                
                
                MM = FastRead(fn{ii});
                indd = find(mean(MM)>mean(MM(:))/8);
                
                GM = nanmean(MM(:,indd),2);
                GM = [];
                obj.Params.QA.out.meanBold{ii} = GM;
                GM1 = diff(GM);
                ind1 = find(abs(GM1)>(std(GM1)*obj.Params.QA.in.GlobalSigThresh))+1;
                
                %%%
                dat = load(obj.Params.Realign.out.realigpars{ii});
                mv = mean(sqrt(sum(diff(dat(:,1:3)).^2,2)));
                sd = std(sqrt(sum(diff(dat(:,1:3)).^2,2)));
                obj.Params.QA.out.meanMV(ii) = mv;
                obj.Params.QA.out.sdMV(ii) = sd;
                
                %%%
                tm = DistMat(dat(:,1:3),0);
                tm = max(tm(:));
                tr = DistMat(dat(:,4:6),0);
                tr = max(tr(:));
                
                %%%
                tmp = diff(dat).^2;
                pos = sqrt(sum(tmp(:,1:3),2));
                ori = sqrt(sum(tmp(:,4:6),2)) * (360/(2*pi));
                ind2 = find(pos>obj.Params.QA.in.MoveThresh)+1;
                ind3 = find(ori>obj.Params.QA.in.RotThresh)+1;
                
                ind = unique([ind1'; ind2; ind3]);
                obj.Params.QA.out.badVols{ii} = ind;
                
                warning off
                obj.Params.QA.out.reasons.SNR(ii,1) = snr(3)<obj.Params.QA.in.SNRthresh;
                obj.Params.QA.out.reasons.meanMV(ii,1) = mv>obj.Params.QA.in.MeanMovementThresh;
                obj.Params.QA.out.reasons.badVols(ii,1) = numel(ind)>obj.Params.QA.in.BadVolThresh;
                obj.Params.QA.out.reasons.totMV(ii,1) = tm>obj.Params.QA.in.TotalMovementThresh;
                obj.Params.QA.out.reasons.totRT(ii,1) = tr>obj.Params.QA.in.TotalRotationThresh;
                warning on;
            end
            obj.Params.QA.out.badruns = find(any(table2array(obj.Params.QA.out.reasons)')');
        end
        
        function obj = proc_SeedMaps(obj,fn,seed_spec,DropVols,movefiles)
            %%% Example:
            %%%     seed_spec = {{'PCC' [  0 -53  26] [10]} ...
            %%%                 {'mPFC' [  0  52 -06] [10]} ...
            %%%                 {'lLPC' [-48 -62  36] [8]}  ...
            %%%                 {'LRPC' [ 46 -62  32] [8]}};
            if ~iscell(fn)
                fn = cellstr(fn);
            end
            
            if nargin<4
                DropVols = false;
            end
            
            if nargin<5
                movefiles = '';
            end
            
            [a bb cc] = fileparts(fn{1});
            outpath = obj.getPath(a,movefiles);
                        
            h = spm_vol([fn{1} ',1']);
            MM = []; vind = [];
            for ii = 1:length(fn)
                M = FastRead(fn{ii});
                if DropVols == 1
                    vind{ii} = setdiff(1:size(M,1),obj.Params.QA.out.badVols{ii});
                    MM = [MM; M(vind{ii},:)];
                else
                    MM = [MM; M];
                    vind{ii} = 1:size(M,1);
                end
            end
            
            ind = find(~isnan(sum(MM)));
            
            %%%
            x = [];
            for zz = 1:numel(seed_spec)
                x{zz} = [];
                
                if numel(seed_spec{zz})==2 && ischar(seed_spec{zz}{2})
                    mask = resizeVol2(spm_vol(seed_spec{zz}{2}),h,[0 0]);
                    i1 = find(mask>0);
                    i1 = intersect(i1,ind);
                    x{zz} = [x{zz}; mean(MM(:,i1),2)];
                    lab{zz} = 'MaskImg.nii';
                elseif numel(seed_spec{zz})==2 && isnumeric(seed_spec{zz}{2})
                    x{zz} = [x{zz}; seed_spec{zz}{2}];
                    lab{zz} = 'PreSpecTS.nii';
                elseif numel(seed_spec{zz}) == 3
                    [ml vi] = getMatCoord(h(1),seed_spec{zz}{2},seed_spec{zz}{3});
                    
                    i1 = intersect(ind,vi);
                    x{zz} = [x{zz}; mean(MM(:,i1),2)];
                    
                    t1 = num2str(seed_spec{zz}{2}); t1 = regexprep([t1], '  ', ' '); t1 = regexprep([t1], '  ', ' '); t1 = regexprep([t1], '  ', ' ');
                    lab{zz} = [num2str(seed_spec{zz}{3}) 'mm_' regexprep([num2str(t1)], ' ', '_') '.nii'];
                end
            end
            
            %%%
            V = h(1);
            vol = nan(V.dim);
            
            for zz = 1:numel(seed_spec)
                b = pinv(x{zz})*MM;
                vol(ind) = atanh(b);
                V.fname = [outpath 'zmap_' seed_spec{zz}{1} '_' bb '_' lab{zz}];
                disp(V.fname);
                V.dt = [16 0];
                spm_write_vol(V,vol);
            end
        end   
         
        function obj = proc_First_Level_Model(obj)
            % Need to run some tests with either missing runs and/or bad
            % runs to make sure the behaviors are appropriate.
            fprintf('\n\n%s\n', 'RUNNING FIRST LEVEL MODEL:');
            wasRun = false;
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            od = pwd;
            % Think about how to utilize the RunIDs field if specified.
            fn = obj.Params.FirsLevelMod.in.fn;
            
            
            [a b c] = fileparts(fn{1});
            if isempty(a); a = pwd; end;
            outpath = obj.getPath(a,[obj.Params.FirsLevelMod.in.movespmfiles filesep obj.Params.FirsLevelMod.in.modelname]);
            
            if ~exist([outpath 'ResMS.nii']);
                wasRun = true;
                
                nscans = [];
                for ii = 1:numel(fn);
                    nscans(ii) = numel(spm_vol(fn{ii}));
                end
                
                runs = 1:numel(fn);
                goodruns = setdiff(runs,obj.Params.QA.out.badruns);
                
                if numel(goodruns) < obj.Params.FirsLevelMod.in.minRuns
                    warning('Not enough good runs to process data.  Stopping now...');
                    return
                end
                
                h = spm_vol([obj.rawfiles{1} ',1']);
                if isempty(obj.Params.FirsLevelMod.in.mp.MicrotimeRes)
                    obj.Params.FirsLevelMod.in.mp.MicrotimeRes = h.dim(3);             %% Changed this from 16 to 30 after fMRI course in Abq. new value = number of slices
                end
                if isempty(obj.Params.FirsLevelMod.in.mp.MicrotimeOnset)
                    obj.Params.FirsLevelMod.in.mp.MicrotimeOnset = floor(h.dim(3)/2);  %% Changed this from 1 to 15 after fMRI course in Abq. new value = middle slice;
                end
                
                h = spm_vol([fn{1} ',1']);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%
                SPM = [];
                SPM.xY.RT = obj.TR;
                SPM.xY.P  = char(fn); %% Will need to subset this where appropriate
                
                %%% Put together SPM.xBF
                spm('defaults','FMRI')
                % The following lines have the side effect of modifying the global
                % defaults variable. This is necessary to pass job.timing.fmri_t to
                % spm_hrf.m. The original values are saved here and restored at the end
                % of this function, after the design has been specified. The original
                % values may not be restored if this function crashes.
                % global defaults
                % olddefs.stats.fmri.fmri_t=spm_get_defaults('stats.fmri.fmri_t');
                % olddefs.stats.fmri.fmri_t0=spm_get_defaults('stats.fmri.fmri_t0');
                % defaults.stats.fmri.t =  P.cm.MicrotimeRes;
                % defaults.stats.fmri.t0 = P.cm.MicrotimeOnset;
                SPM.xBF.UNITS =          obj.Params.FirsLevelMod.in.mp.Units;
                SPM.xBF.dt =             obj.TR/obj.Params.FirsLevelMod.in.mp.MicrotimeRes;
                SPM.xBF.T     =          obj.Params.FirsLevelMod.in.mp.MicrotimeRes;
                SPM.xBF.T0    =          obj.Params.FirsLevelMod.in.mp.MicrotimeOnset;
                
                if strcmpi('hrf',obj.Params.FirsLevelMod.in.mp.basis_func)
                    SPM.xBF.name = 'hrf';
                elseif strcmpi('hrf (with time derivative)',obj.Params.FirsLevelMod.in.mp.basis_func)
                    SPM.xBF.name = 'hrf (with time derivative)';
                elseif strcmpi('fir',obj.Params.FirsLevelMod.in.mp.basis_func)
                    SPM.xBF.name = 'Finite Impulse Response';
                    SPM.xBF.length  = obj.Params.FirsLevelMod.in.mp.length;
                    SPM.xBF.order   = obj.Params.FirsLevelMod.in.mp.order;
                    %keyboard
                end
                
                
                
                SPM.xBF = spm_get_bf(SPM.xBF);
                
                SPM.xBF.Volterra = obj.Params.FirsLevelMod.in.mp.Volterra;
                
                %keyboard;
                SPM.nscan = sum(nscans);
                
                SPM.xGX.iGXcalc = obj.Params.FirsLevelMod.in.mp.iGXcalc;
                SPM.xGX.sGXcalc = obj.Params.FirsLevelMod.in.mp.sGXcalc;
                SPM.xGX.sGMsca =  obj.Params.FirsLevelMod.in.mp.sGMsca;
                
                SPM.xVi.form = obj.Params.FirsLevelMod.in.mp.AR;
                
                for ii = 1:numel(goodruns)
                    SPM.xX.K(ii).HParam = obj.Params.FirsLevelMod.in.mp.HP_filt;
                end
                
                %%%
                dat = obj.RunInfo;
                SPM.Sess = [];
                c = 0;
                for ii = goodruns;
                    c = c+1;
                    SPM.nscan(c) = nscans(ii);
                    i1 = find(dat.Run == ii);
                    
                    uni = unique(dat.Condition,'stable');
                    
                    if strcmpi(obj.Params.FirsLevelMod.in.mp.Units ,'scans')
                        mx = nscans(ii);
                    else
                        mx = nscans(ii)*obj.TR;
                    end
                    
                    cc = 0;
                    for jj = 1:length(uni);
                        i2 = contains(uni{jj}, dat.Condition(i1));
                        if isempty(i2);
                            continue;
                        end
                        
                        if all(dat.Onset(i1(i2))<mx)
                            cc = cc+1;
                            SPM.Sess(c).U(cc).name{1} = uni{jj};
                            SPM.Sess(c).U(cc).ons =     dat.Onset(i1(i2));
                            SPM.Sess(c).U(cc).dur =     dat.Duration(i1(i2));
                            SPM.Sess(c).U(cc).P.name =  obj.Params.FirsLevelMod.in.mp.TempMod;
                            SPM.Sess(c).U(cc).P.h = 0;      %%?
                            SPM.Sess(c).C.C = [];           %%?
                            SPM.Sess(c).C.name = cell(0);   %%?
                        else
                            error('Onsets go beyond the end of the run');
                        end
                    end
                end
                %%%
                
                c = 0;
                if obj.Params.FirsLevelMod.in.mp.addBadVolRegs == 1;
                    disp('Scrubbing Volumes');
                    if isfield(obj.Params.QA.out,'badVols')
                        for ii = goodruns
                            c = c+1;
                            tmp = obj.Params.QA.out.badVols{ii};
                            X = zeros(nscans(ii),numel(tmp));
                            for jj = 1:numel(tmp);
                                X(tmp(jj),jj) = 1;
                                SPM.Sess(c).C.name{end+1} = ['bv_' num2str(jj)];
                            end
                            SPM.Sess(c).C.C = X;
                        end
                    end
                end
                
                if obj.Params.FirsLevelMod.in.mp.addMotionRegressors == 1;
                    disp('Adding motion parameters to the model');
                    c = 0;
                    for ii = goodruns
                        c = c+1;
                        mot = load(obj.Params.Realign.out.realigpars{ii});
                        for kk = 1:size(mot,2);
                            SPM.Sess(c).C.name{end+1} = ['mot_' num2str(kk)];
                        end
                        SPM.Sess(c).C.C = [SPM.Sess(c).C.C mot];
                    end
                end
                cd(outpath);
                
                save('SPM.mat','SPM');
                
                delete('beta*');
                delete('mask*');
                delete('ResMS*');
                delete('RPV*');
                delete('SPM.mat');
                delete('ess*');
                delete('con*');
                delete('spmF*');
                delete('spmT*');
                
%                 try
                SPM = spm_fmri_spm_ui(SPM); %% Will cause compiled script to fail.
%                 catch; keyboard; end
                
                if obj.Params.FirsLevelMod.in.mp.DisableThresholdMasking == 1;
                    SPM.xM.TH(:) = -Inf;  %% disable threshold masking
                else
                    SPM.xM.I = 1;
                end
                
                if isfield(obj.Params.FirsLevelMod.in.mp, 'ExplicitMask');
                    SPM.xM.VM = obj.Params.FirsLevelMod.in.mp.ExplicitMask;
                end
                SPM = spm_spm(SPM);
                
            else
                disp('First Level Model Estimation is Complete.');
            end
            
            if obj.Params.FirsLevelMod.in.mp.ScreenTaskCorrMot && (~isfield(obj.Params.FirsLevelMod.out,'tcm') || isempty(obj.Params.FirsLevelMod.out.tcm));
                if ~exist('SPM','var')
                    load([outpath 'SPM.mat']);
                end
                goodruns = setdiff(1:numel(obj.rawfiles),obj.Params.QA.out.badruns);
                nscans = SPM.nscan;
                
                dat = obj.RunInfo;
                disp('Screening for task correlated motion');
                conds = unique(dat.Condition);
                ccc = 0; r2 = [];
                tmpscans = nscans;
                r2 = nan(numel(obj.rawfiles),numel(conds));
                for ii = goodruns;
                    ccc = ccc+1;
                    rows =  1+sum(tmpscans(1:ccc-1)):sum(tmpscans(1:ccc));
                    mot = load(obj.Params.Realign.out.realigpars{ii});
                    
                    xx = [];
                    for jj = 1:length(conds);
                        th1 = contains(['.*' num2str(ccc) '.*' conds{jj} '.*bf\(1\).*'],SPM.xX.name);
                        if isempty(th1);
                            xx(:,jj) = zeros(numel(rows),1);
                            continue;
                        else
                            xx(:,jj) = SPM.xX.X(rows,th1);
                        end
                    end
                    
                    pred = demean(mot)*(pinv(demean(mot))*demean(xx));
                    r2(ii,:) = diag(corr(pred,xx).^2);
                end
                disp('% Variance explained for each condition for each run');
                disp([conds'; num2cell(r2)])
                %obj.Params.FirsLevelMod.out.tcm = [conds'; num2cell(r2)];
                obj.Params.FirsLevelMod.out.tcm = array2table(r2,'VariableNames',conds');
                % What should be done with this information?  Could try to
                % use this to modify the bad runs list, and then call the
                % method again to restart the analysis
            end
            
            cd(od);
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            spmfile = [outpath 'SPM.mat'];
            obj.Params.FirsLevelMod.out.path = outpath;
            obj.Params.FirsLevelMod.out.spm = spmfile;
            
            obj.UpdateHist(obj.Params.FirsLevelMod,'proc_First_Level_Model',spmfile,wasRun);
            fprintf('\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function obj = proc_Gen_Contrasts(obj)
            % And then there's this guy!  Will also need to update/edit the
            % CreateConVec Script.
            fprintf('\n%s\n', 'Creating Contrast Images:');
            wasRun = false;
             
            od = pwd;
            outpath = [fileparts(obj.Params.FirsLevelMod.out.spm) filesep];
            groupcons = obj.Params.Contrasts.in.GroupConFold;
            
            if 0>inf;%numel(dir_wfp([outpath 'spmT*.nii']))>0
                load(obj.Params.FirsLevelMod.out.spm);
                obj.Params.Contrasts.out.con = [];
                for ii = 1:numel(SPM.xCon)
                    obj.Params.Contrasts.out.con.(SPM.xCon(ii).name) = [outpath SPM.xCon(ii).Vcon.fname];
                end
                disp('Contrasts have already been generated');
            else
                cd(outpath);
                wasRun = true;

                delete('ess*');
                delete('con*');
                delete('spmF*');
                delete('spmT*');
                
                x = load('SPM.mat');
                SPM = x.SPM;
                Con = [];
                cc = 0;
                        
                cons = obj.Params.Contrasts.in.con;
                c = 0;
                for ii = 1:length(cons);
                    [vec L R] = CreateConVec2(cons(ii).left,cons(ii).right,SPM,cons(ii).WeightWithin,cons(ii).BlockThresh);
                    if all(isnan(vec)); continue; end
                    if isempty(vec)
                        warning('how does this happen?')
                        keyboard;
                        % continue
                    end
                    
                    obj.Params.Contrasts.out.conInfo(ii).con = cons(ii);
                    obj.Params.Contrasts.out.conInfo(ii).L = L;
                    obj.Params.Contrasts.out.conInfo(ii).R = R;
                    obj.Params.Contrasts.out.conInfo(ii).vec = vec;
                    obj.Params.Contrasts.out.conInfo(ii).note = '';
                    
                    if (~isempty(L.counts) && sum(L.counts)<obj.Params.Contrasts.in.MinEvents) || (~isempty(R.counts) && sum(R.counts)<obj.Params.Contrasts.in.MinEvents)
                        obj.Params.Contrasts.out.conInfo(ii).note = ['Skipping contrast "' cons(ii).name '": too few events.'];
                        disp(obj.Params.Contrasts.out.conInfo(ii).note);
                        continue
                    end
                    
                    cc = cc+1;
                    Con(cc).name = obj.Params.Contrasts.in.con(ii).name;
                    if size(vec,1)==1
                        Con(cc).STAT = 'T';
                    else
                        Con(cc).STAT = 'F';
                    end
                    Con(cc).c = vec';
                end
                                
                if length(Con)==0
                    error('No Contrasts met the Minimum Event# criteria');
                end
                
                
                for ii = 1:length(Con)
                    xCon(ii) = spm_FcUtil('Set',Con(ii).name,Con(ii).STAT,'c',Con(ii).c,SPM.xX.xKXs);
                end
                
                SPM.xCon = xCon;
                SPM = spm_contrasts(SPM,1:length(SPM.xCon));
                
                load SPM.mat
                for ii = 1:length(SPM.xCon)
                    curname = SPM.xCon(ii).Vcon.fname(1:end-4);
                    newname = SPM.xCon(ii).name;
                    SPM.xCon(ii).Vcon.fname = [curname(1:4) newname '.nii'];
                    
                    movefile([curname '.nii'], [curname(1:4) newname '.nii']);
                    
                    if SPM.xCon(ii).STAT == 'T';
                        copyfile(['spmT_' sprintf('%0.4i',ii) '.nii'], ['spmT_' newname '.nii']);
                    end
                    if SPM.xCon(ii).STAT == 'F';
                        copyfile(['spmF_' sprintf('%0.4i',ii) '.nii'], ['spmF_' newname '.nii']);
                    end
                    
                    try
                    obj.Params.Contrasts.out.con.(SPM.xCon(ii).name) = [outpath SPM.xCon(ii).Vcon.fname];
                    catch
                        keyboard;
                    end
                end
                save SPM.mat SPM
                
                if obj.Params.Contrasts.in.MoveContrasts == 1
                    disp('Copying Contrast Images');
                    
                    
                    %%% Remove any existing contrasts for this session;
                    p = groupcons;
                    p = [':' p];
                    ind = find(p==':');
                    
                    for ii = 1:length(ind)-1;
                        tmp = p(ind(ii)+1:ind(ii+1)-1);
                        tmp = [tmp filesep fn '*'];
                        try; ls(tmp); end
                        delete(tmp);
                    end
                    
                    
                    for ii = 1:length(SPM.xCon)
                        if exist([groupcons SPM.xCon(ii).name])==0
                            mkdir([groupcons SPM.xCon(ii).name]);
                        end
                        
                        copyfile(SPM.xCon(ii).Vcon.fname, [groupcons SPM.xCon(ii).name filesep obj.sessionname '_' SPM.xCon(ii).name '.nii']);
                    end
                end
            end
            
            if isempty(obj.Params.Contrasts.out.conInfo);
                cons = obj.Params.Contrasts.in.con;
                for ii = 1:length(cons);
                    [vec L R] = CreateConVec2(cons(ii).left,cons(ii).right,SPM,cons(ii).WeightWithin,cons(ii).BlockThresh);
                    
                    obj.Params.Contrasts.out.conInfo(ii).con = cons(ii);
                    obj.Params.Contrasts.out.conInfo(ii).L = L;
                    obj.Params.Contrasts.out.conInfo(ii).R = R;
                    obj.Params.Contrasts.out.conInfo(ii).vec = vec;
                    obj.Params.Contrasts.out.conInfo(ii).note = '';
                    
                    if (~isempty(L.counts) && sum(L.counts)<obj.Params.Contrasts.in.MinEvents) || (~isempty(R.counts) && sum(R.counts)<obj.Params.Contrasts.in.MinEvents)
                        obj.Params.Contrasts.out.conInfo(ii).note = ['Skipping contrast "' cons(ii).name '": too few events.'];
                        disp(obj.Params.Contrasts.out.conInfo(ii).note);
                        continue
                    end
                end
            end
            
            cd(od);
            flds = fields(obj.Params.Contrasts.out.con);
            fn = obj.Params.Contrasts.out.con.(flds{1});
            obj.UpdateHist(obj.Params.FirsLevelMod,'proc_Gen_Contrasts',fn,wasRun);
            fprintf('\n');
        end
        
        function obj = extractTS(obj,fn, seeds, BV)
            %%%  Extract a time series from a 4D Nifti file.
            %%%
            %%% INPUTS:
            %%% fn:     Path to 4D nifti file.
            %%% seeds:  Cell Array of Inputs { {'Label1' [mni coords] [seed diameter]} ...}
            %%% opt:    If opt = 1, the time series will be extracted.
            %%%         If opt = 0, only the matrix indcies of the seed will be
            %%%         computed.
            %%%
            %%% OUTPUTS
            %%% D.fn =  A new filename based on the location and size of the seed.
            %%% D.vox_loc = location of the seed in 3D matrix coordinates.
            %%% D.vec_loc = location of the seed in vector coordinates.
            %%% D.ts =  The time series as computed by the first eigen vector method.
            %%% D.tsm = The time series as computed by a simple mean of the values in
            %%% the seed at each time point.
            %%%
            
            obj.Params.Seeds.out.seeds = [];
            
            data = table;
           
            for ii = 1:numel(fn)
                badVols = obj.Params.QA.out.badVols{ii};
                
                V1 = spm_vol(fn{ii});
                
                [tmp XYZmm] = spm_read_vols(V1(1));
                
                ind = setdiff(1:numel(V1),badVols);
                M2 = FastRead(V1(ind));
                
                for zz = 1:length(seeds)
                    obj.Params.Seeds.out.seeds{ii} = seeds{ii};
                    
                    if isnumeric(seeds{zz}{2});
                        rad = seeds{zz}{3}/2;
                        dist = sqrt(sum((XYZmm' - repmat(seeds{zz}{2},size(XYZmm,2),1)).^2,2));
                        ind = find(dist <= rad);
                        if numel(ind)==0;
                            ind = find(dist==min(dist));
                        end
                    elseif ischar(seeds{zz}{2});
                        h = spm_vol(seeds{zz}{2});
                        mi = resizeVol(h,V1(1));
                        ind = find(mi==1);
                        [aa1 aa2] = fileparts(seeds{zz}{2});
                    end
                    
                    
                    dd = XYZmm(:,ind)';
                    clear dd2;
                    [dd2(:,1) dd2(:,2) dd2(:,3)] = ind2sub(V1(1).dim,ind);
                    
                    y = M2(:,ind);
                    data.(seeds{zz}{1}) = nanmean(y,2);
                end
                obj.Params.Seeds.out.data{ii} = data;
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
                info = regexprep(sprintf('%-30s%s', [cmd ':'], ['Last run by ' theUser ' on ' theTime]),'\n','');
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
        
       %replaced with openIMGv2.m ==> 
%         function [ out_fn isgz ] = check_gz(obj,ext,fn)
%             if strcmp(ext,'.gz')
%                 fprintf(['\n Gunzipping: ' fn ]);
%                 out_fn=cell2char(gunzip(fn));
%                 fprintf('...done!\n');
%                 isgz=true;
%             else
%                 out_fn=fn{ii};
%                 isgz=false;
%             end
%         end
        
       %replaced with writeIMG.m  ==> 
%         function [ out_fn ] = do_gz(fn)
%             fprintf(['\n Gzipping: ' fn ' ...']);
%             gzip(fn);
%             out_fn=[fn '.gz' ] ;
%             fprintf ('...done\n');
%             system(['rm ' fn])
%         end
        
    end
end


