classdef dwiMRI_Session  < dynamicprops & matlab.mixin.SetGet
    %%% Written by:
    %%%             Rodrigo Perea (rpereacamargo@mgh.harvard.edu)
    %%%             Aaron Schultz (aschultz@martinos.org)
    %%
    %%      Dependencies (and tested in):
    %%          -FreeSurfer v6.0
    %%          -SPM8
    %%          -Ants tools 1.7.1
    %%          -DSI_studio_vApr19_2017
    %%          -FSL 5.0.9
    %%
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
        exec_onecmd = '' ; %Sets the opt option when initialize to execute one single instruction.
        sessionname = '';
        rawfiles = '';
        root=''; %to be modified if specified as an argument.
        dosave = false;
        wasLoaded = false;
        rawstructural = '';
        %%% Change things so that these are the only version used.
        fsdir = '';
        fsubj = '';
        FSdata='';
        objectHome = '';
        T1 = '' ;
        projectID = ''
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
        Trkland = [] ;
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
            catch
                disp('Cant get collectiondate')
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
            obj.Params.DCM2NII.specific_vols = [];
            obj.Params.DCM2NII.in.nvols = [];
            obj.Params.DCM2NII.in.fsl2std_matfile = '';
            obj.Params.DCM2NII.in.fsl2std_param = '';
            
            obj.Params.DCM2NII.out.location = [];
            obj.Params.DCM2NII.out.fn = [];
            obj.Params.DCM2NII.out.bvecs='';
            obj.Params.DCM2NII.out.bvals='';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.DropVols.in.tmin = ''; %char type as they will be called in system('')
            obj.Params.DropVols.in.tsize = '';
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
            obj.Params.GradNonlinCorrect.in.fn = '';
            obj.Params.GradNonlinCorrect.in.fslroi= [0 1];
            obj.Params.GradNonlinCorrect.out.b0='';
            obj.Params.GradNonlinCorrect.out.warpfile = [];
            obj.Params.GradNonlinCorrect.out.meannii = [];
            obj.Params.GradNonlinCorrect.out.fn = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.B0MoCo.FS = '';
            obj.Params.B0MoCo.in.movefiles = '';
            obj.Params.B0MoCo.in.prefix = 'moco_';
            obj.Params.B0MoCo.in.nDoF = '12' ;
            obj.Params.B0MoCo.in.sh_rotate_bvecs = '';
            
            obj.Params.B0MoCo.in.fn = '';
            obj.Params.B0MoCo.in.bvals = '';
            obj.Params.B0MoCo.in.bvecs = '';
            
            obj.Params.B0MoCo.out.fn = '';
            obj.Params.B0MoCo.out.bvals = '';
            obj.Params.B0MoCo.out.bvecs = '';
            
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
            obj.Params.EddyMotion.in.fn_eddy = '';
            obj.Params.EddyMotion.in.fn_motion = '';
            
            obj.Params.EddyMotion.out.fn_motion = '';
            obj.Params.EddyMotion.out.vals.initb0_mean = [];
            obj.Params.EddyMotion.out.vals.initb0_std = [];
            obj.Params.EddyMotion.out.vals.rel_mean = [];
            obj.Params.EddyMotion.out.vals.rel_std = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.B0mean.in.movefiles = '';
            obj.Params.B0mean.in.fn='';
            obj.Params.B0mean.in.b0_nvols=[];
            obj.Params.B0mean.out.fn='';
            obj.Params.B0mean.out.allb0s='';
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.MaskAfterEddy.in.movefiles = '';
            obj.Params.MaskAfterEddy.in.fn = '';
            obj.Params.MaskAfterEddy.in.prefix = '';
            obj.Params.MaskAfterEddy.in.b0 = '';
            obj.Params.MaskAfterEddy.out.initmask = '';
            obj.Params.MaskAfterEddy.out.finalmask = ''; %Corrected for inconsistent brain edges value on dtifit after using the option --wls
            obj.Params.MaskAfterEddy.out.brainonly = '';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.CoRegMultiple.in.fn = '' ;
            obj.Params.CoRegMultiple.in.b0 = '' ;
            obj.Params.CoRegMultiple.in.bvals = '' ;
            obj.Params.CoRegMultiple.in.bvecs = '' ;
            obj.Params.CoRegMultiple.in.movefiles = '' ;
            obj.Params.CoRegMultiple.in.ref_iteration = [] ; % All images will be registered to this iteration (in ADRC, 7p5_set1, index 1 is for 2p7_set4!)
            obj.Params.CoRegMultiple.in.ref_prefix = '' ;
            obj.Params.CoRegMultiple.in.ref_file = '' ;
            
            obj.Params.CoRegMultiple.out.fn = '' ;
            obj.Params.CoRegMultiple.out.bvals = '' ;
            obj.Params.CoRegMultiple.out.bvecs = '' ;
            obj.Params.CoRegMultiple.out.matfile = '' ;
            obj.Params.CoRegMultiple.out.combined_fn = '' ;
            obj.Params.CoRegMultiple.out.combined_bvals = '' ;
            obj.Params.CoRegMultiple.out.combined_bvecs = '' ;
            obj.Params.CoRegMultiple.out.combined_b0 = '' ;
            obj.Params.CoRegMultiple.out.combined_bet = '' ;
            obj.Params.CoRegMultiple.out.combined_mask = '' ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.Dtifit.in.movefiles = '';
            obj.Params.Dtifit.in.fn='';
            obj.Params.Dtifit.in.bvecs='';
            obj.Params.Dtifit.in.bvals='';
            obj.Params.Dtifit.in.mask='';
            obj.Params.Dtifit.in.prefix = '' ;
            obj.Params.Dtifit.out.prefix='';
            obj.Params.Dtifit.out.FA='';
            obj.Params.Dtifit.out.RD = '';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.GQI.in.movefiles = '';
            obj.Params.GQI.in.prefix = '' ;
            obj.Params.GQI.in.fn = '' ;
            obj.Params.GQI.in.bvecs = '';
            obj.Params.GQI.in.bvals = '';
            obj.Params.GQI.in.mask = '' ;
            
            obj.Params.GQI.in.method = '4'; %for gqi
            obj.Params.GQI.in.num_fiber = '3';  %modeling 3 fiber population
            obj.Params.GQI.in.param0 = '1.25';
            
            
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
            obj.Params.AntsReg.in.threads = '1' ;
            obj.Params.AntsReg.in.transform = 's' ;
            obj.Params.AntsReg.in.radius = '4' ;
            obj.Params.AntsReg.in.precision = 'd' ;
            obj.Params.AntsReg.in.prefix = '';
            obj.Params.AntsReg.out.fn = '';
            obj.Params.AntsReg.out.FA = '';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.FROIS2dwi.in.dim = '3';
            obj.Params.FROIS2dwi.in.threads = '1';
            obj.Params.FROIS2dwi.in.transform = 's';
            obj.Params.FROIS2dwi.in.radius = '4';
            obj.Params.FROIS2dwi.in.precision = 'd' ;
            obj.Params.FROIS2dwi.in.ref  = '';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            obj.Params.Skeletonize.out.diffmetrics={ 'FA' 'RD' 'AxD' 'MD' } ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.Skel_TOI.in.location = '' ;
            obj.Params.Skel_TOI.in.masks = '' ;
            obj.Params.Skel_TOI.in.ext = '.nii.gz' ;
            obj.Params.Skel_TOI.out = [] ; %this should be populated with all the masked TOIs
            obj.Params.Skel_TOI.in.suffix = '';
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.FreeSurfer.dir = '' ;
            obj.Params.FreeSurfer.in.T1 = '' ;
            obj.Params.FreeSurfer.in.T2 = '' ;
            obj.Params.FreeSurfer.in.T2exist = false ;
            
            obj.Params.FreeSurfer.out.aparcaseg = '' ;
            obj.Params.FreeSurfer.init_location = '';
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.FS2dwi.in.movefiles = '';
            obj.Params.FS2dwi.in.b0 = '' ;
            obj.Params.FS2dwi.in.aparcaseg = '' ;
            obj.Params.FS2dwi.in.tmpfile_aparcaseg = '' ;
            obj.Params.FS2dwi.in.aparcaseg2009 = '' ;
            obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = '' ;
            obj.Params.FS2dwi.in.hippofield_left = '' ;
            obj.Params.FS2dwi.in.hippofield_right = '' ;
            obj.Params.FS2dwi.in.tmpfile_hippo = '' ;
            
            obj.Params.FS2dwi.out.xfm_dwi2FS = '' ;
            obj.Params.FS2dwi.out.fn_aparc = '';
            obj.Params.FS2dwi.out.fn_aparc2009 = '';
            
            %             %%%%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%%%
            %             %%ALL THESE BELOW ARE INSTANCES OF A PREVIOUS CLASS (USED FOR
            %             %%CODE RECYLING...)
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%% BEGIN Things %%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
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
        function obj = proc_dcm2nii(obj,~,dest,out_filename)
            fprintf('\n\n%s\n', ['CONVERTING DCM TO NII (Total:  ' num2str(numel(obj.Params.DCM2NII.in)) ' volumes)']);
            %%% MODULE IS COMPLETE
            wasRun = false;
            nfn = [];
            for ii=1:numel(obj.Params.DCM2NII.in)
                if nargin>2 && ~isempty(dest)
                    obj.Params.DCM2NII.out(ii).location = dest;
                end
                
                if nargin>3 && ~isempty(out_filename)
                    obj.Params.DCM2NII.out(ii).out_filename = out_filename;
                end
                
                %Shortening naming comventions:
                clear in_file out_file
                in_file=strtrim([obj.dcm_location filesep obj.Params.DCM2NII.in(ii).first_dcmfiles]);
                out_file=obj.Params.DCM2NII.out(ii).fn;
                outpath=obj.Params.DCM2NII.out(ii).location;
                obj.Params.DCM2NII.in(ii).fsl2std_matfile = [outpath 'fsl2std.matfile'];
                
                obj.Params.DCM2NII.out(ii).bvecs=strrep(obj.Params.DCM2NII.out(ii).fn,'.nii.gz','.voxel_space.bvecs');
                obj.Params.DCM2NII.out(ii).bvals=strrep(obj.Params.DCM2NII.out(ii).fn,'.nii.gz','.bvals');
                %Create out_file directory if doesnt exist:
                if ~exist(obj.Params.DCM2NII.out(ii).location,'dir')
                    clear exec_cmd
                    exec_cmd = ['mkdir -p ' obj.Params.DCM2NII.out(ii).location ];
                    obj.RunBash(exec_cmd);
                end
                
                %Processing starts here:
                if exist(in_file,'file') ~= 0 %check if in_file exists
                    if exist(out_file,'file') == 0 %check if out_file exists
                        %Check whether we get the specific number of volumes:
                        if  obj.Params.DCM2NII.specific_vols == obj.Params.DCM2NII.in(ii).nvols
                            clear exec_cmd
                            exec_cmd=['mri_convert ' in_file ' ' out_file ];
                            obj.RunBash(exec_cmd,44);
                            
                            fprintf('\nFslreorienting to standard...')
                            %Reorient to std that nii files -->
                            exec_cmd=['fslreorient2std ' out_file ' ' out_file ];
                            obj.RunBash(exec_cmd);
                            %%%
                            
                            %Reorient to std the bvecs -->
                            if isempty(obj.Params.DCM2NII.in(ii).fsl2std_param) %due to inconsitencies with bvecs from mri_convert, this should be initialize in the child class
                                exec_cmd=['fslreorient2std ' out_file ' > ' obj.Params.DCM2NII.in(ii).fsl2std_matfile ];
                                obj.RunBash(exec_cmd);
                            else
                                exec_cmd=['echo -e ''' obj.Params.DCM2NII.in(ii).fsl2std_param ''' > ' obj.Params.DCM2NII.in(ii).fsl2std_matfile ];
                                obj.RunBash(exec_cmd);
                                
                            end
                            %%%
                            
                            %Now dealing with bvecs:
                            disp('Fslreorienting the bvecs now...')
                            temp_bvec=[outpath 'temp.bvec' ];
                            exec_cmd=[obj.init_rotate_bvecs_sh ' ' ...
                                ' ' obj.Params.DCM2NII.out(ii).bvecs ...
                                ' ' obj.Params.DCM2NII.in(ii).fsl2std_matfile ...
                                ' ' temp_bvec  ];
                            obj.RunBash(exec_cmd);
                            
                            system(['mv ' temp_bvec ' ' obj.Params.DCM2NII.out(ii).bvecs ]);
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
        
        function obj = proc_drop_vols(obj)
            fprintf('\n%s\n', 'PERFORMING PROC_DROP_VOLS():');
            wasRun=false;
            for ii=1:numel(obj.Params.DropVols.in.fn)
                clear cur_fn;
                if iscell(obj.Params.DropVols.in.fn{ii})
                    cur_fn=cell2char(obj.Params.DropVols.in.fn{ii});
                else
                    cur_fn=obj.Params.DropVols.in.fn{ii};
                end
                [a b c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.DropVols.in.movefiles);
                clear outfile
                obj.Params.DropVols.out.fn{ii} = [ outpath obj.Params.DropVols.in.prefix  b c ];
                
                %Droppping volumes in the DWIs (niftii):
                if exist( obj.Params.DropVols.out.fn{ii},'file')==0
                    fprintf(['\n Dropping volumes  (fslroi <input> <output> ' ...
                        obj.Params.DropVols.in.tmin ' ' obj.Params.DropVols.in.tsize  ') ...' ...
                        'Iteration: ' num2str(ii) ]);
                    exec_cmd = [ 'fslroi ' obj.Params.DropVols.in.fn{ii} ...
                        ' ' obj.Params.DropVols.out.fn{ii}  ' ' obj.Params.DropVols.in.tmin ...
                        ' ' obj.Params.DropVols.in.tsize ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n ');
                    wasRun=true;
                    obj.UpdateHist(obj.Params.DropVols,'proc_drop_vols()', obj.Params.DropVols.out.fn{ii},wasRun);
                else
                    fprintf(['Dropping vols for: ' b c ' is complete. \n']) ;
                end
                %Droppping volumes in the DWIs (bvals):
                obj.Params.DropVols.in.bvals{ii}=strrep(obj.Params.DropVols.in.fn{ii},'.nii.gz','.bvals');
                obj.Params.DropVols.out.bvals{ii}=strrep(obj.Params.DropVols.out.fn{ii},'.nii.gz','.bvals');
                if exist( obj.Params.DropVols.out.bvals{ii},'file')==0
                    exec_cmd=(['sed 1,' obj.Params.DropVols.in.tmin 'd ' ...
                        obj.Params.DropVols.in.bvals{ii} ' > ' obj.Params.DropVols.out.bvals{ii} ]);
                    obj.RunBash(exec_cmd);
                end
                %Droppping volumes in the DWIs (bvecs):
                obj.Params.DropVols.in.bvecs{ii}=strrep(obj.Params.DropVols.in.fn{ii},'.nii.gz','.voxel_space.bvecs');
                obj.Params.DropVols.out.bvecs{ii}=strrep(obj.Params.DropVols.out.fn{ii},'.nii.gz','.voxel_space.bvecs');
                if exist( obj.Params.DropVols.out.bvecs{ii},'file')==0
                    exec_cmd=(['sed 1,' obj.Params.DropVols.in.tmin 'd ' ...
                        obj.Params.DropVols.in.bvecs{ii} ' > ' obj.Params.DropVols.out.bvecs{ii} ]);
                    obj.RunBash(exec_cmd);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%OPTIONAL:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Check if outfiles are in row format. If so, put them in col
            %format
            [a , b ] = size(obj.Params.DropVols.out.bvals);
            if a < b
                obj.Params.DropVols.out.bvals = obj.Params.DropVols.out.bvals';
            end
            [a , b ] = size(obj.Params.DropVols.out.bvecs);
            if a < b
                obj.Params.DropVols.out.bvecs = obj.Params.DropVols.out.bvecs';
            end
            %%%%%%%%%%%%%%%%%%END OF OPTIONAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
        end
        
        function obj = proc_gradient_nonlin_correct(obj)
            fprintf('\n%s\n', 'PERFORMING PROC_GRADIENT_NONLIN_CORRECT():');
            wasRun = false;
            in_fn = obj.Params.GradNonlinCorrect.in.fn;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %APPLYING GRADIENT_NON_LIN CORRECTION TO set4 and THEN EXPAND
            %IT TO ALL OTHER IMAGES
            %this will take set4 in ADRC data (shouldn;t affect it. Need to double-check!)
            [a b c ] = fileparts(in_fn{1});
            outpath=obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);
            
            %first extract the first b0s (if it doesn't exist):
            obj.Params.GradNonlinCorrect.in.b0{1}=[outpath 'firstb0_' b c ];
            if exist(obj.Params.GradNonlinCorrect.in.b0{1},'file')==0
                fprintf(['\nExtracting the first b0 of: ' obj.Params.GradNonlinCorrect.in.b0{1} ]);
                exec_cmd=['fslroi ' in_fn{1} ' ' obj.Params.GradNonlinCorrect.in.b0{1} ...
                    ' ' num2str(obj.Params.GradNonlinCorrect.in.fslroi) ];
                obj.RunBash(exec_cmd);
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
                exec_cmd=['sh ' obj.sh_gradfile ' '  first_b0_infile ' ' first_b0_outfile ' ' gradfile ' '];
                obj.RunBash(exec_cmd,44);
                wasRun = true;
            else
                [~, bb, cc ] = fileparts(first_b0_outfile);
                fprintf([ 'The gnc file ' bb cc  ' is complete.\n'])
            end
            
            %%Apply the correction to the first_b0
            if exist(first_b0_outfile,'file')==0
                exec_cmd=['applywarp -i ' first_b0_infile ' -r ' first_b0_infile ...
                    ' -o ' first_b0_outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile{1} ...
                    ' --interp=spline' ];
                fprintf(['\nGNC: Applying warp to first_b0_file: '  first_b0_infile]);
                obj.RunBash(exec_cmd);
                
                fprintf('...done\n');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %NOW APPLYING IT TO ALL IMAGES:
            %%% Apply the correction to all the subsequent diffusion images.
            for ii=1:numel(in_fn)
                [a b c ] = fileparts(in_fn{ii});
                obj.Params.GradNonlinCorrect.out.fn{ii}=[ outpath 'gnc_' b c ];
                if exist(obj.Params.GradNonlinCorrect.out.fn{ii},'file')==0
                    fprintf('\nGNC: Applying warp field to all the other images:\n');
                    fprintf(['~~~~> ' obj.Params.GradNonlinCorrect.out.fn{ii} '.']);
                    exec_cmd = ['applywarp -i ' obj.Params.GradNonlinCorrect.in.fn{ii} ' -r ' first_b0_infile ...
                        ' -o ' obj.Params.GradNonlinCorrect.out.fn{ii} ' -w ' obj.Params.GradNonlinCorrect.out.warpfile{1} ' --interp=spline'];
                    obj.RunBash(exec_cmd);
                    fprintf('....done\n');
                    wasRun = true;
                    obj.UpdateHist(obj.Params.GradNonlinCorrect,'proc_gradient_nonlin_correct',obj.Params.GradNonlinCorrect.out.fn{ii},wasRun);
                end
            end
            
            %%%%%%%%%%%%%%%%  END OF OPTIONAL: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Check if outfiles are in row format. If so, put them in col
            %format
            [a , b ] = size(obj.Params.GradNonlinCorrect.out.fn);
            if a < b
                obj.Params.GradNonlinCorrect.out.fn = obj.Params.GradNonlinCorrect.out.fn';
            end
            %%%%%%%%%%%%%%%%  END OF OPTIONAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function obj = proc_bet2(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_BET2():');
            
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
                    exec_cmd = [ 'bet2 ' obj.Params.Bet2.in.fn{ii} ' ' obj.Params.Bet2.out.skull{ii}  ' -m -f ' num2str(obj.Params.Bet2.in.fracthrsh) ];
                    obj.RunBash(exec_cmd);
                    system(['mv ' obj.Params.Bet2.out.skull{ii} '_mask.nii.gz ' obj.Params.Bet2.out.mask{ii} ] ) ;
                    wasRun=true;
                    obj.UpdateHist(obj.Params.Bet2,'proc_bet2', obj.Params.Bet2.out.mask{ii},wasRun);
                else
                    [aa, bb, cc ] = fileparts(obj.Params.Bet2.out.mask{ii});
                    fprintf([' File ' bb cc ' is now comple. \n']) ;
                end
            end
        end
        
        function obj = proc_eddy(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_EDDY():');
            
            for ii=1:numel(obj.Params.Eddy.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Eddy.in.fn{ii})
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
                    obj.RunBash(exec_cmd);
                    fprintf(' ...done\n');
                end
                %(dependency) Attempting to create index file:
                obj.Params.Eddy.out.fn_index{ii}= [ outpath 'index.txt' ] ;
                if exist(obj.Params.Eddy.out.fn_index{ii},'file')==0
                    fprintf(['\n Creating ' obj.Params.Eddy.out.fn_index{ii}  ]);
                    exec_cmd=['echo " ' num2str(obj.Params.Eddy.in.index) ' " >> ' obj.Params.Eddy.out.fn_index{ii}  ];
                    obj.RunBash(exec_cmd);
                    fprintf(' ...done\n');
                end
                %Attempting to run eddy_openmp now:
                obj.Params.Eddy.out.fn{ii} = [ outpath obj.Params.Eddy.in.prefix  b c ];
                obj.Params.Eddy.out.bvecs{ii} = [ outpath obj.Params.Eddy.in.prefix strrep(b,'.nii','.eddy_rotated_bvecs') ];
                if exist( obj.Params.Eddy.out.fn{ii},'file')==0
                    try
                        fprintf(['\nApplying eddy in: ' obj.Params.Eddy.in.fn{ii} ]);
                        fprintf('\n this will take a couple of minutes...');
                        exec_cmd=[ 'eddy_openmp --imain=' obj.Params.Eddy.in.fn{ii} ...
                            ' --mask=' obj.Params.Eddy.in.mask{ii} ...
                            ' --index=' obj.Params.Eddy.out.fn_index{ii} ...
                            ' --acqp='  obj.Params.Eddy.out.fn_acqp{ii}  ...
                            ' --bvecs='  obj.Params.Eddy.in.bvecs{ii} ...
                            ' --bvals=' obj.Params.Eddy.in.bvals{ii}  ...
                            ' --repol --out=' [ outpath obj.Params.Eddy.in.prefix strrep(b,'.nii','') ]  ];
                        obj.RunBash(exec_cmd,44);
                        fprintf('...done \n');
                        
                        wasRun=true;
                        obj.UpdateHist(obj.Params.Eddy,'proc_eddy', obj.Params.Eddy.out.fn{ii},wasRun);
                    catch
                        errormsg=['PROC_EDDY: Cannnot run eddy in: ' ...
                            obj.Params.Eddy.in.fn{ii} 'please double check parameters?\n' ];
                        obj.UpdateErrors(errormsg);
                    end
                else
                    [aa, bb, cc] = fileparts(obj.Params.Eddy.out.fn{ii});
                    fprintf(['File ' bb cc ' is now complete \n']) ;
                end
            end
        end
        
        function obj = proc_b0s_MoCo(obj)
            %Motion Correction for interspersed b0s.
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING MOTION CORRECTION - PROC_B0S_MOCO():');
            %Check if the FreeSurfer location exists (due to bbreg dependency):
            if exist( obj.Params.B0MoCo.FS,'dir') ~=0 %if so continue, else break!
                for jj=1:numel(obj.Params.B0MoCo.in.fn)
                    clear cur_fn;
                    %Variable type fixing issues:
                    if iscell(obj.Params.B0MoCo.in.fn{jj})
                        cur_fn=cell2char(obj.Params.B0MoCo.in.fn{jj});
                    else
                        cur_fn=obj.Params.B0MoCo.in.fn{jj};
                    end
                    %Splitting the naming convention:
                    [a b c ] = fileparts(cur_fn);
                    outpath=obj.getPath(a,obj.Params.B0MoCo.in.movefiles);
                    
                    
                    %Initializing obj.Params.B0MoCo.out.XX:
                    obj.Params.B0MoCo.out.fn{jj}= [ outpath obj.Params.B0MoCo.in.prefix b c ] ;
                    obj.Params.B0MoCo.out.bvecs{jj}= [ outpath  obj.Params.B0MoCo.in.prefix strrep(b,'.nii','.voxel_space.bvecs') ];
                    obj.Params.B0MoCo.out.bvals{jj}= [ outpath  obj.Params.B0MoCo.in.prefix strrep(b,'.nii','.bvals') ];
                    
                    %Splitting the current DWI:
                    clear tmp_fslsplit;
                    tmp_fslsplit=[ outpath 'tmp' filesep ...
                        'tmp_' strrep(b,'.nii','') filesep ];
                    if exist([tmp_fslsplit '0000.nii.gz' ],'file') == 0
                        exec_cmd = (['mkdir -p ' tmp_fslsplit ] );
                        obj.RunBash(exec_cmd);
                        exec_cmd=(['fslsplit ' obj.Params.B0MoCo.in.fn{jj} ...
                            ' ' tmp_fslsplit ' -t ']);
                        obj.RunBash(exec_cmd);
                    end
                    
                    %Record bval information:
                    clear tmp_bval_idx
                    tmp_bval_idx=load(obj.Params.B0MoCo.in.bvals{jj});
                    
                    
                    
                    %%%%%%%%STEP 1: EXRTACT INFO FROM B0s%%%%%%%%%%%%%%%%%%
                    flag_b0_idx = [] ; %this will denote what b0idx to use (later in dwi_b0_idx{ii} variable (at end of for loop)
                    
                    
                    if exist(obj.Params.B0MoCo.out.fn{jj}, 'file') == 0
                        %Apply bbreg to all b0s:
                        for ii=1:numel(tmp_bval_idx)
                            dwi_idx=ii-1  ; %indexing from fslsplit starts at 0 not 1 so we need this additional indexing
                            %Get single dwi (indexed at 0000 or 0010)
                            if ii < 11 % (idx start at 1 in matlab but fslsplit idx starts at 0, hence the difference!!
                                cur_in_dwi{ii} = [ tmp_fslsplit '000' num2str(dwi_idx) '.nii.gz' ];
                            else
                                cur_in_dwi{ii} = [ tmp_fslsplit '00' num2str(dwi_idx) '.nii.gz' ];
                            end
                            [ ~, bn_curdwi{ii}, ~ ] = fileparts(cur_in_dwi{ii});
                            
                            
                            
                            %INIT VARIABLE (ISOLATED from other inits (that are
                            %inside the if b0 loop because we will use them to
                            %apply warp in step 2 of this method...
                            out_dwi2firstb0{ii} =  [ tmp_fslsplit 'dwi' strrep(bn_curdwi{ii},'.nii','') '_2_modfirstb0.nii.gz' ];
                            
                            %Bbregister starts here (only for b0s!):
                            if tmp_bval_idx(ii) == 0 %it it's a b0 image
                                flag_b0_idx = ii ;
                                %Init variables (b0 dependable):
                                out_trnii{ii}= [ tmp_fslsplit 'dwi_' strrep(bn_curdwi{ii},'.nii','') '_2_fsT1.nii.gz' ];
                                out_dwi2fslmat{ii} = [ tmp_fslsplit 'mat_dwi' strrep(bn_curdwi{ii},'.nii','') '_2_fsT1.mat' ];
                                out_dwi2firstmat{ii} = [ tmp_fslsplit 'mat_dwi' strrep(bn_curdwi{ii},'.nii','') '_2_firstdwi.mat' ];
                                out_dwi2firstmat_rot_only{ii} = [ tmp_fslsplit 'rotonly_mat_dwi' strrep(bn_curdwi{ii},'.nii','') '_2_firstdwi.mat' ];
                                out_dwi2fsmat{ii} =  [ tmp_fslsplit 'fs_mat_dwi' strrep(bn_curdwi{ii},'.nii','') '_2_fsT1.dat' ];
                                
                                
                                
                                
                                %bbregister on every b0:
                                if exist(out_dwi2fslmat{ii},'file') == 0
                                    exec_cmd = ([ 'export SUBJECTS_DIR=' obj.FS_location ';' ...
                                        ' bbregister --s ' obj.sessionname ' --' obj.Params.B0MoCo.in.nDoF ...
                                        ' --mov ' cur_in_dwi{ii} ' --dti --o ' out_trnii{ii}  ...
                                        ' --init-header --reg ' out_dwi2fsmat{ii} ...
                                        ' --fslmat ' out_dwi2fslmat{ii}   ]);
                                    obj.RunBash(exec_cmd,44);
                                end
                                %Convert dwi2T1 to T12dwi and apply_warp
                                %(different from the 1st b0 and consequent
                                %ones, hence the if else statement:
                                if ii == 1
                                    out_fsl2firstdwimat = [ tmp_fslsplit 'mat_fslT1_2_dwi' strrep(bn_curdwi{1},'.nii','') '.mat' ];
                                    %Convert dwi2T1 into T12dwi:
                                    if exist(out_fsl2firstdwimat,'file') == 0
                                        exec_cmd=[ 'convert_xfm -omat ' out_fsl2firstdwimat ...
                                            ' -inverse ' out_dwi2fslmat{ii} ];
                                        obj.RunBash(exec_cmd);
                                    end
                                    %Apply the disco_rel here using applywarp
                                    if exist(out_dwi2firstb0{ii} ,'file') == 0
                                        exec_cmd=(['applywarp -i ' cur_in_dwi{ii} ' -r ' cur_in_dwi{1} ...
                                            ' -o ' out_dwi2firstb0{ii}  ' -w ' obj.Params.B0MoCo.in.grad_rel{1} ...
                                            ' --rel --interp=spline' ]);
                                        obj.RunBash(exec_cmd);
                                    end
                                    
                                    %Creating the eye.mat function for
                                    %reference on the firstb0 for rotation
                                    %matrices that will be used for rot. bvecs
                                    if exist(out_dwi2firstmat_rot_only{ii}, 'file' ) == 0
                                        exec_cmd=(['echo -e "1 0 0 0 \n0 1 0 0 \n0 0 1 0 \n0 0 0 1 " > ' ...
                                            out_dwi2firstmat_rot_only{ii} ]);
                                        obj.RunBash(exec_cmd);
                                    end
                                    
                                else
                                    %Convert dwi2T1 into T12dwi (concating the firstb0):
                                    if exist(out_dwi2firstmat{ii},'file') == 0
                                        exec_cmd=[ 'convert_xfm -omat ' out_dwi2firstmat{ii}  ...
                                            ' -concat ' out_fsl2firstdwimat ' ' out_dwi2fslmat{ii} ];
                                        obj.RunBash(exec_cmd);
                                    end
                                    
                                    %Apply the disco_rel here using applywarp
                                    if exist(out_dwi2firstb0{ii} ,'file') == 0
                                        exec_cmd=(['applywarp -i ' cur_in_dwi{ii} ' -r ' cur_in_dwi{1} ...
                                            ' -o ' out_dwi2firstb0{ii} ' -w ' obj.Params.B0MoCo.in.grad_rel{1} ...
                                            ' --postmat=' out_dwi2firstmat{ii} ' --interp=spline --rel ' ]);
                                        obj.RunBash(exec_cmd);
                                    end
                                    
                                    %Extracing the rotation matrix only using
                                    %avscale (this will be used for modifying
                                    %the bvecs output)
                                    if exist(out_dwi2firstmat_rot_only{ii}, 'file' ) == 0
                                        exec_cmd=(['avscale ' out_dwi2firstmat{ii} ...
                                            ' | head -5 | tail -4 > ' out_dwi2firstmat_rot_only{ii} ]);
                                        obj.RunBash(exec_cmd);
                                    end
                                end
                            end
                            dwi_b0_idx{ii} = flag_b0_idx ; %This will keep the idx of what b0 to be used (for when applying warp to all b0s)
                        end
                        
                        %%%%%%%%STEP 2: APPLY INFO DERIVED FROM B0s%%%%%%%%%%%%
                        %Now applying all the information derived from B0s to
                        %all the DWIs (including the B0s)
                        
                        disp('In proc_b0s_MoCo Step 2: apply info from derived b0s');
                        disp(['In file iteration number: ' num2str(jj) ]);
                        fprintf('\n In DWI: ')
                        for ii=1:numel(tmp_bval_idx)
                            %Applying warp to all non-B0s (as it is done in
                            %Step1)
                            
                            if exist(out_dwi2firstb0{ii},'file') == 0
                                fprintf([ '...' bn_curdwi{ii} ] );
                                exec_cmd=(['applywarp -i ' cur_in_dwi{ii} ...
                                    ' -r ' cur_in_dwi{1} ' -o ' out_dwi2firstb0{ii} ...
                                    ' -w ' obj.Params.B0MoCo.in.grad_rel{1} ...
                                    ' --postmat=' out_dwi2firstmat{dwi_b0_idx{ii}} ...
                                    ' --interp=spline' ]);
                                obj.RunBash(exec_cmd);
                            end
                            
                            %Modify bvecs information now...
                            AA=1;
                        end
                        fprintf('...done');
                        
                        
                        
                        %%%%%%%%STEP 3: APPLY ROT_MATFILES TO BVECS%%%%%%%%
                        if exist(obj.Params.B0MoCo.out.bvecs{jj},'file') == 0
                            fprintf('Applying rotation only .mat files to bvecs..');
                            TEMP_BVEC = load(obj.Params.B0MoCo.in.bvecs{jj});
                            for pp=1:size(TEMP_BVEC,1)
                                exec_cmd=[obj.Params.B0MoCo.in.sh_rotate_bvecs ...
                                    ' ' num2str(TEMP_BVEC(pp,:)) ...
                                    ' ' out_dwi2firstmat_rot_only{dwi_b0_idx{pp}} ...
                                    ' >> ' (obj.Params.B0MoCo.out.bvecs{jj}) ];
                                obj.RunBash(exec_cmd);
                            end
                            fprintf('...done');
                        end
                        
                        %Copy bval_in to bvals_out (just a simple copy as
                        %its not affected):
                        if exist(obj.Params.B0MoCo.out.bvals{jj},'file') == 0
                            exec_cmd=['cp ' obj.Params.B0MoCo.in.bvals{jj} ' ' obj.Params.B0MoCo.out.bvals{jj} ];
                            obj.RunBash(exec_cmd);
                        end
                        
                        
                        %%%%%%%CREATING OUTPUT DWI NOW %%%%%%%%%%%%%%%%%%%%
                        %Fslmerging all the corrected values at this point
                        all_DWI_str = [];
                        for kk=1:numel(out_dwi2firstb0)
                            %First check that all volumes exist:
                            if ~exist(out_dwi2firstb0{kk},'file')==2
                                error_txt=[ 'Error: ' out_dwi2firstb0{kk} ' does not exist!. Please check' ];
                                exec_cmd=['echo "' error_txt '"  >> ' outpath 'error.txt' ]
                                obj.RunBash(exec_cmd);
                            end
                            all_DWI_str = [ all_DWI_str ' ' out_dwi2firstb0{kk} ]  ;
                        end
                        exec_cmd = [ 'fslmerge -t ' obj.Params.B0MoCo.out.fn{jj}  ' ' all_DWI_str ];
                        fprintf('\nfslmerging all newer b0MoCo corrected volumes...')
                        obj.RunBash(exec_cmd);
                        wasRun=true;
                        obj.UpdateHist(obj.Params.B0MoCo,'proc_B0MoCo', obj.Params.B0MoCo.out.fn{jj},wasRun);
                        
                        fprintf('...done \n')
                        
                    else
                        [aa, bb, cc ] = fileparts(obj.Params.B0MoCo.out.fn{jj});
                        fprintf([ 'File ' bb cc ' is complete.\n' ]);
                    end
                    
                    
                end
            end
        end
        
        function obj = proc_get_eddymotion(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING MOTION EXTRACTION (from eddy) - PROC_GET_EDDYMOTION():');
            
            for ii=1:numel(obj.Params.EddyMotion.in.fn_eddy)
                clear cur_fn;
                %Dealing with INPUT:
                %obj.Params.EddyMotion.in.fn_motion --> using restricted
                %movement file instead of the non_restricted as it give us
                %the actual movement disregarding the tranlation in the PE
                %direction (double check:
                %https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide )
                obj.Params.EddyMotion.in.fn_motion{ii} = ...
                    strrep(obj.Params.EddyMotion.in.fn_eddy{ii},'.nii.gz','.eddy_restricted_movement_rms');
                if exist(obj.Params.EddyMotion.in.fn_motion{ii} )
                    if iscell(obj.Params.EddyMotion.in.fn_motion{ii})
                        cur_fn=cell2char(obj.Params.EddyMotion.in.fn_motion{ii});
                    else
                        cur_fn=obj.Params.EddyMotion.in.fn_motion{ii};
                    end
                    [a b c ] = fileparts(cur_fn);
                else
                    error(['No eddy_motion file exist with this name: ' ...
                        obj.Params.EddyMotion.in.fn_motion{ii} ' exiting...']);
                end
                %Dealing with OUTPUT:
                %Attempting to create the eddy_motion file:
                outpath=obj.getPath(a,obj.Params.EddyMotion.in.movefiles);
                obj.Params.EddyMotion.out.fn_motion{ii}= [ outpath 'motion_values.txt' ] ;
                if exist(obj.Params.EddyMotion.out.fn_motion{ii},'file')==0
                    fprintf(['\nGetting motion parameters from: ' obj.Params.EddyMotion.in.fn_motion{ii} ]);
                    %Break the command so I can denote the newline (\n
                    %character):
                    cmd_1=sprintf([ 'awk ''{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} \n ' ...
                        '  END {for (i=1;i<=NF;i++) { \n' ...
                        ' printf "' ] );
                    cmd_2= ('%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)} ' ) ;
                    cmd_3= sprintf(['\n }'' ' obj.Params.EddyMotion.in.fn_motion{ii}  ' > ' obj.Params.EddyMotion.out.fn_motion{ii} ]);
                    exec_cmd = [ cmd_1 cmd_2 cmd_3 ] ;
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                    wasRun=true;
                else
                    [~, bb, cc ] =  fileparts(obj.Params.EddyMotion.out.fn_motion{ii});
                    fprintf(['File ' bb cc 'exists\n']) ;
                end
                %Populating the variables needed:
                if isempty(obj.Params.EddyMotion.out.vals.initb0_mean)
                    [ ~, obj.Params.EddyMotion.out.vals.initb0_mean{ii} ]  = ...
                        system([ 'cat ' obj.Params.EddyMotion.out.fn_motion{ii} ' | head -1 | awk '' {print $1} '' ' ]);
                end
                if isempty(obj.Params.EddyMotion.out.vals.initb0_std)
                    [ ~, obj.Params.EddyMotion.out.vals.initb0_std{ii} ]  = ...
                        system([ 'cat ' obj.Params.EddyMotion.out.fn_motion{ii} ' | head -1 | awk '' {print $2} '' ' ]);
                end
                if isempty(obj.Params.EddyMotion.out.vals.rel_mean)
                    [ ~ , obj.Params.EddyMotion.out.vals.rel_mean{ii} ] = ...
                        system([ 'cat ' obj.Params.EddyMotion.out.fn_motion{ii} ' | tail -1 | awk '' {print $1} '' ' ]);
                end
                if isempty( obj.Params.EddyMotion.out.vals.rel_std)
                    [ ~ , obj.Params.EddyMotion.out.vals.rel_std{ii} ] = ...
                        system([ 'cat ' obj.Params.EddyMotion.out.fn_motion{ii} ' | tail -1 | awk '' {print $2} '' ' ]);
                end
                
                obj.UpdateHist(obj.Params.Eddy,'proc_eddy_motion', obj.Params.Eddy.out.fn{ii},wasRun);
            end
        end
        
        function obj = proc_meanb0(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_MEANB0():');
            
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
                        obj.RunBash(exec_cmd);
                        fprintf('...done');
                    end
                    if exist( obj.Params.B0mean.out.fn{ii},'file')==0
                        fprintf(['\n Meaning all b0s from : ' cur_fn]);
                        exec_cmd=[ 'fslmaths ' obj.Params.B0mean.out.allb0s{ii} ...
                            ' -Tmean ' obj.Params.B0mean.out.fn{ii}];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.B0mean,'proc_b0mean', obj.Params.B0mean.out.fn{ii},wasRun);
                    else
                        [~, bb, cc] = fileparts(obj.Params.B0mean.out.fn{ii});
                        fprintf(['The file ' bb cc ' is complete. \n']);
                    end
                catch
                    errormsg=['PROC_B0MEAN: Cannnot create the following meanB0 from:'  ...
                        cur_fn 'Please check this input location!\n' ];
                    obj.UpdateErrors(errormsg);
                    
                end
            end
        end
        
        function obj = proc_mask_after_eddy(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_MASK_AFTER_EDDY():');
            
            for ii=1:numel(obj.Params.MaskAfterEddy.in.fn)
                [a b c ] = fileparts(obj.Params.MaskAfterEddy.in.fn{ii});
                outpath=obj.getPath(a,obj.Params.MaskAfterEddy.in.movefiles);
                clear outfile
                
                %Check whether the input was a b0 or the whole DWI:
                [~ , tmp_nvols ] =system(['fslinfo ' obj.Params.MaskAfterEddy.in.fn{ii} ' | grep ^dim4 | awk ''{print $2}''' ]);
                if strcmp(strtrim(tmp_nvols),'1')
                    %Here, we assume you pass a b0 argument (no need for
                    %more complicated naming) --> compatible with dwi_HAB.m
                    obj.Params.MaskAfterEddy.in.b0{ii} =  obj.Params.MaskAfterEddy.in.fn{ii} ;
                    obj.Params.MaskAfterEddy.out.brainonly{ii} = [ outpath obj.Params.MaskAfterEddy.in.prefix  '_brainonly.nii.gz' ];
                else
                    %Else we need extra steps to extract the first b0...
                    %(compatible with dwi_ADRC.m)
                    obj.Params.MaskAfterEddy.in.b0{ii} = [ outpath 'b0_' b c ] ;
                    exec_cmd = (['fslroi ' obj.Params.MaskAfterEddy.in.fn{ii} ...
                        ' '  obj.Params.MaskAfterEddy.in.b0{ii} ' 0 1 ' ]);
                    if exist(obj.Params.MaskAfterEddy.in.b0{ii},'file') == 0
                        fprintf(['\n proc_mask_after_eddy: Extracting the first b0 for iteration: ' num2str(ii) ]);
                        obj.RunBash(exec_cmd);
                        fprintf('...done\n');
                    end
                    obj.Params.MaskAfterEddy.out.brainonly{ii} = strrep(obj.Params.MaskAfterEddy.in.b0{ii},'b0_','brainonly_');
                end
                
                obj.Params.MaskAfterEddy.out.initmask{ii} = strrep(obj.Params.MaskAfterEddy.out.brainonly{ii},'brainonly_','initmask_');
                obj.Params.MaskAfterEddy.out.finalmask{ii} = strrep(obj.Params.MaskAfterEddy.out.brainonly{ii},'brainonly_','finalmask_');
                if exist( obj.Params.MaskAfterEddy.out.finalmask{ii},'file')==0
                    fprintf(['\nExtracting the brain only using bet2 for : ' obj.Params.MaskAfterEddy.in.b0{ii} ]);
                    %Initial mask creted:
                    exec_cmd = [ 'bet2 ' obj.Params.MaskAfterEddy.in.b0{ii} ' ' obj.Params.MaskAfterEddy.out.brainonly{ii}  ' -m -f ' num2str(obj.Params.MaskAfterEddy.in.fracthrsh) ];
                    obj.RunBash(exec_cmd);
                    exec_cmd=(['mv ' obj.Params.MaskAfterEddy.out.brainonly{ii} '_mask.nii.gz ' obj.Params.MaskAfterEddy.out.initmask{ii} ] ) ;
                    obj.RunBash(exec_cmd);
                    %Final mask created (for --wls in dtifit!)
                    exec_cmd = [ 'fslmaths ' obj.Params.MaskAfterEddy.in.b0{ii} ' -thr 10 -bin -mas ' ...
                        obj.Params.MaskAfterEddy.out.initmask{ii} ' '  obj.Params.MaskAfterEddy.out.finalmask{ii} ];
                    obj.RunBash(exec_cmd);
                    wasRun=true;
                    obj.UpdateHist(obj.Params.Bet2,'proc_mask_after_eddy', obj.Params.MaskAfterEddy.out.finalmask{ii},wasRun);
                else
                    [aa, bb, cc] = fileparts(obj.Params.MaskAfterEddy.out.finalmask{ii});
                    fprintf(['File ' bb cc ' is now complete. \n']) ;
                end
            end
        end
        
        %Use when multiple DWIs sequences are acquired
        function obj = proc_coreg_multiple(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_COREG_MULTIPLE():');
            %First, select the ref_files needed:
            for ii=1:numel(obj.Params.CoRegMultiple.in.fn)
                [a b c ] = fileparts(obj.Params.CoRegMultiple.in.fn{ii});
                outpath=obj.getPath(a,obj.Params.CoRegMultiple.in.movefiles);
                %Copy bvecs, bvals and fn
                obj.Params.CoRegMultiple.out.fn{ii} = [ outpath 'coreg_' b c ];
                obj.Params.CoRegMultiple.out.bvals{ii} = [ outpath 'coreg_' strrep(b,'.nii','.bvals') ] ;
                obj.Params.CoRegMultiple.out.bvecs{ii} = [ outpath 'coreg_' strrep(b,'.nii','.bvecs') ] ;
                
                if ii ==  obj.Params.CoRegMultiple.in.ref_iteration
                    tmp_b_split = strsplit(b,'_');
                    obj.Params.CoRegMultiple.in.ref_prefix = cell2char(strrep(tmp_b_split(end),'.nii',''));
                    obj.Params.CoRegMultiple.in.ref_file = obj.Params.CoRegMultiple.in.b0{ii};
                    %Copy the files in this for loop (since nothing will be done to ref)
                    if exist(obj.Params.CoRegMultiple.out.fn{ii},'file') == 0
                        %*.nii.gz:
                        exec_cmd=(['cp ' obj.Params.CoRegMultiple.in.fn{ii} ...
                            ' ' obj.Params.CoRegMultiple.out.fn{ii}  ]);
                        obj.RunBash(exec_cmd);
                    end
                    if exist(obj.Params.CoRegMultiple.out.bvals{ii},'file') == 0
                        %*.bvals:
                        exec_cmd=(['cp ' obj.Params.CoRegMultiple.in.bvals{ii} ...
                            ' ' obj.Params.CoRegMultiple.out.bvals{ii}  ]);
                        obj.RunBash(exec_cmd);
                    end
                    if exist(obj.Params.CoRegMultiple.out.bvecs{ii},'file') == 0
                        %*.bvecs:
                        exec_cmd=(['sh ' obj.col2rows_sh ' ' obj.Params.CoRegMultiple.in.bvecs{ii} ...
                            ' > ' obj.Params.CoRegMultiple.out.bvecs{ii}  ]);
                        obj.RunBash(exec_cmd);
                    end
                end
            end
            %Now apply flirt in all the other images:
            for ii=1:numel(obj.Params.CoRegMultiple.in.fn)
                %Initializing matfiles:
                [aa bb cc ] = fileparts(obj.Params.CoRegMultiple.in.fn{ii});
                tmp_bb_split = strsplit(bb,'_');
                obj.Params.CoRegMultiple.out.matfile{ii} = [ outpath ...
                    'xfm_b0' cell2char(strrep(tmp_bb_split(end),'.nii','_2_')) ...
                    obj.Params.CoRegMultiple.in.ref_prefix '.mat' ];
                %End of init
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ii ~=  obj.Params.CoRegMultiple.in.ref_iteration
                    %Dealing with *.nii.gz:
                    
                    if exist(obj.Params.CoRegMultiple.out.matfile{ii},'file') == 0
                        %Creating matfiles:
                        exec_cmd=['flirt -in ' obj.Params.CoRegMultiple.in.b0{ii} ...
                            ' -ref ' obj.Params.CoRegMultiple.in.ref_file ...
                            ' -omat ' obj.Params.CoRegMultiple.out.matfile{ii} ...
                            ' -interp spline ' ];
                        fprintf(['\n proc_coreg_multiple: rigid coregistratio with iteration: ' num2str(ii)]);
                        obj.RunBash(exec_cmd);
                        fprintf('...done\n');
                    end
                    if exist(obj.Params.CoRegMultiple.out.fn{ii},'file') == 0
                        %Applying the matfile from b0s:
                        exec_cmd=(['applywarp -i ' obj.Params.CoRegMultiple.in.fn{ii} ' -r ' obj.Params.CoRegMultiple.in.ref_file ...
                            ' -o ' obj.Params.CoRegMultiple.out.fn{ii} ...
                            ' --postmat=' obj.Params.CoRegMultiple.out.matfile{ii} ' --interp=spline ' ]);
                        fprintf(['\n proc_coreg_multiple: applying warp (iter ' num2str(ii) ')']);
                        obj.RunBash(exec_cmd);
                        fprintf('...done\n');
                    end
                    
                    %Dealing with *.bvecs:
                    for tohide=1:1
                        if exist(obj.Params.CoRegMultiple.out.bvecs{ii},'file') == 0
                            %*.bvecs:
                            %% from rows 3-by-XX to cols XX-by-3:
                            out_tmp_bvecs{ii}= [outpath 'tmp_bvecs_iter' num2str(ii) ] ;
                            exec_cmd=(['sh ' obj.col2rows_sh ' ' obj.Params.CoRegMultiple.in.bvecs{ii} ...
                                ' > ' out_tmp_bvecs{ii}  ]);
                            obj.RunBash(exec_cmd);
                            
                            
                            %Extracing the rotation matrix only using
                            %avscale (this will be used for modifying
                            %the bvecs output)
                            tmp_rot_avscale{ii}= [outpath 'tmp_rotonly_iter' num2str(ii) ] ;
                            if exist(tmp_rot_avscale{ii}, 'file' ) == 0
                                exec_cmd=(['avscale ' obj.Params.CoRegMultiple.out.matfile{ii} ...
                                    ' | head -5 | tail -4 > ' tmp_rot_avscale{ii} ]);
                                obj.RunBash(exec_cmd);
                            end
                            
                            
                            %Now apply rotation values to the newer bvecs:
                            fprintf('Applying rotation only .mat files to bvecs..');
                            TEMP_BVEC{ii} = load(out_tmp_bvecs{ii});
                            
                            %remove cause it will be replaced:
                            if exist(obj.Params.CoRegMultiple.out.bvecs{ii},'file') == 2
                                system(['rm ' obj.Params.CoRegMultiple.out.bvecs{ii} ]);
                            end
                            for pp=1:size(TEMP_BVEC{ii},1)
                                exec_cmd=[obj.b0MoCo_rotate_bvecs_sh ...
                                    ' ' num2str(TEMP_BVEC{ii}(pp,:)) ...
                                    ' ' tmp_rot_avscale{ii}  ...
                                    ' >> ' (obj.Params.CoRegMultiple.out.bvecs{ii}) ];
                                obj.RunBash(exec_cmd);
                            end
                            system(['rm ' tmp_rot_avscale{ii}]);
                            system(['rm ' out_tmp_bvecs{ii}]);
                            fprintf('...done');
                            
                        end
                    end
                    %Copying bvals:
                    if exist(obj.Params.CoRegMultiple.out.bvals{ii},'file') == 0
                        %*.bvals:
                        exec_cmd=(['cp ' obj.Params.CoRegMultiple.in.bvals{ii} ...
                            ' ' obj.Params.CoRegMultiple.out.bvals{ii}  ]);
                        obj.RunBash(exec_cmd);
                    end
                    
                end
            end
            
            %Creating the combined niftiis:
            for tohide=1:1
                obj.Params.CoRegMultiple.out.combined_fn = [outpath 'combined_preproc_' num2str(numel(obj.Params.CoRegMultiple.out.fn)) 'sets' '.nii.gz' ] ;
                if exist( obj.Params.CoRegMultiple.out.combined_fn , 'file') == 0
                    tmp_nii_cmd = [] ;
                    for ii=1:numel(obj.Params.CoRegMultiple.out.fn)
                        if ii ~=  obj.Params.CoRegMultiple.in.ref_iteration
                            tmp_nii_cmd = [  tmp_nii_cmd ' ' obj.Params.CoRegMultiple.out.fn{ii} ] ;
                        end
                    end
                    exec_cmd = [ 'fslmerge -t ' obj.Params.CoRegMultiple.out.combined_fn ...
                        ' ' obj.Params.CoRegMultiple.out.fn{obj.Params.CoRegMultiple.in.ref_iteration} ...
                        ' ' tmp_nii_cmd ];
                    fprintf('\nMerging niis...')
                    obj.RunBash(exec_cmd);
                    fprintf('...done\n')
                end
                
                %Combining bvals:
                obj.Params.CoRegMultiple.out.combined_bvals =  [outpath 'combined_preproc_' num2str(numel(obj.Params.CoRegMultiple.out.bvals)) 'sets' '.bvals' ] ;
                if exist( obj.Params.CoRegMultiple.out.combined_bvals , 'file') == 0
                    tmp_bvals_cmd = [] ;
                    for ii=1:numel(obj.Params.CoRegMultiple.out.bvals)
                        if ii ~=  obj.Params.CoRegMultiple.in.ref_iteration
                            tmp_bvals_cmd = [  tmp_bvals_cmd ' ' obj.Params.CoRegMultiple.out.bvals{ii} ] ;
                        end
                    end
                    exec_cmd = [ 'cat '  obj.Params.CoRegMultiple.out.bvals{obj.Params.CoRegMultiple.in.ref_iteration} ...
                        ' ' tmp_bvals_cmd ' > ' obj.Params.CoRegMultiple.out.combined_bvals     ];
                    fprintf('\nMerging bvals...')
                    obj.RunBash(exec_cmd);
                    fprintf('\n...done.')
                end
                
                %Combining bvecs:
                obj.Params.CoRegMultiple.out.combined_bvecs =  [outpath 'combined_preproc_' num2str(numel(obj.Params.CoRegMultiple.out.bvecs)) 'sets' '.bvecs' ] ;
                if exist(obj.Params.CoRegMultiple.out.combined_bvecs , 'file') == 0
                    tmp_bvecs_cmd = [] ;
                    for ii=1:numel(obj.Params.CoRegMultiple.out.bvecs)
                        if ii ~=  obj.Params.CoRegMultiple.in.ref_iteration
                            tmp_bvecs_cmd = [  tmp_bvecs_cmd ' ' obj.Params.CoRegMultiple.out.bvecs{ii} ] ;
                        end
                    end
                    exec_cmd = [ 'cat '  obj.Params.CoRegMultiple.out.bvecs{obj.Params.CoRegMultiple.in.ref_iteration} ...
                        ' ' tmp_bvecs_cmd ' > ' obj.Params.CoRegMultiple.out.combined_bvecs     ];
                    fprintf('\nMerging bvecs...')
                    obj.RunBash(exec_cmd);
                    fprintf('\n...done.')
                end
            end
            
            %Extracting the combined b0/mask:
            obj.Params.CoRegMultiple.out.combined_b0 = [outpath 'combined_preproc_b0.nii.gz'] ;
            obj.Params.CoRegMultiple.out.combined_bet = [outpath 'combined_preproc_bet.nii.gz'] ;
            
            obj.Params.CoRegMultiple.out.combined_mask = [outpath 'combined_preproc_bet_mask.nii.gz'] ;
            %b0:
            if exist(obj.Params.CoRegMultiple.out.combined_b0,'file') == 0
                exec_cmd = [ 'cp '  obj.Params.CoRegMultiple.in.b0{obj.Params.CoRegMultiple.in.ref_iteration} ...
                    ' ' obj.Params.CoRegMultiple.out.combined_b0 ]
                fprintf('\n Copying b0 combined...')
                obj.RunBash(exec_cmd);
                fprintf('..done \n')
            end
            %bet and mask:
            if exist(obj.Params.CoRegMultiple.out.combined_bet,'file') == 0
                exec_cmd = [ 'bet2 ' obj.Params.CoRegMultiple.out.combined_b0 ...
                    ' ' obj.Params.CoRegMultiple.out.combined_bet  ' -f 0.3 -m ' ]
                fprintf('\n Bet masking combined...')
                obj.RunBash(exec_cmd);
                fprintf('..done \n')
                movefile([obj.Params.CoRegMultiple.out.combined_bet '_mask.nii.gz'],  obj.Params.CoRegMultiple.out.combined_mask);
            end
            fprintf(' proc_coreg_multiple(obj) is complete\n');
        end
        
        function obj = proc_dtifit(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_DTIFIT():');
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
                obj.Params.Dtifit.out.FA{ii} = [ outpath  obj.Params.Dtifit.in.prefix '_FA.nii.gz' ] ;
                obj.Params.Dtifit.out.prefix{ii} = [ outpath  obj.Params.Dtifit.in.prefix ];
                % try
                %Attempting to dtifit:
                if exist( obj.Params.Dtifit.out.FA{ii},'file')==0
                    
                    %Check if the *.nii exists but not *.nii.gz
                    if exist(strrep(obj.Params.Dtifit.out.FA{ii},'.nii.gz','.nii'))==0
                        fprintf('\nDtifit reconstruction...');
                        exec_cmd=[ 'dtifit -k ' obj.Params.Dtifit.in.fn{ii} ...
                            ' -o ' obj.Params.Dtifit.out.prefix{ii} ...
                            ' -m ' obj.Params.Dtifit.in.mask{ii} ...
                            ' -r ' obj.Params.Dtifit.in.bvecs{ii} ...
                            ' -b ' obj.Params.Dtifit.in.bvals{ii} ...
                            ' --wls --sse' ]; %weighted least squared for improving inadequate noisy data
                        %REF:
                        %
                        obj.RunBash(exec_cmd,44);
                        fprintf('...done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.Eddy,'proc_dtifit', obj.Params.Dtifit.out.FA{ii},wasRun);
                        
                        %Coopy L1 to AxD:
                        exec_cmd = [ 'cp  ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L1') ' ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','AxD') ] ;
                        obj.RunBash(exec_cmd);
                    else
                        display([ 'Gzipping: ' strrep(obj.Params.Dtifit.out.FA{ii},'.nii.gz','.nii')  '...'] );
                        system([ 'gzip ' strrep(obj.Params.Dtifit.out.FA{ii},'.nii.gz','.nii') ] );
                    end
                else
                    [aa,bb,cc] = fileparts(obj.Params.Dtifit.out.FA{ii} );
                    fprintf([ ' File ' bb cc ' is now complete.\n'])
                end
                %                 catch
                %                     errormsg=['PROC_DTIFIT: Cannot apply dtifit:'  ...
                %                         'Please check dtifit input location!\n' ];
                %                     error('Exiting now, cannot apply dtifit');
                %                     obj.UpdateErrors(errormsg);
                %                 end
                
                %Outputting RD:
                obj.Params.Dtifit.out.RD{ii} = strrep(obj.Params.Dtifit.out.FA{ii},'FA','RD');
                if exist(obj.Params.Dtifit.out.RD{ii})==0
                    
                    exec_cmd=[ ' fslmaths ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L2') ...
                        ' -add ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','L3') ...
                        ' -div 2 ' obj.Params.Dtifit.out.RD{ii}  ];
                    fprintf('\nCreating RD dtifit...');
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
        end
        
        function obj = proc_gqi(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_GQI():');
            for ii=1:numel(obj.Params.GQI.in.fn)
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
                obj.Params.GQI.out.btable{ii}= [ outpath  obj.Params.GQI.in.prefix '_btable.txt' ] ;
                
                %Attempting to create b_table:
                if exist(obj.Params.GQI.out.btable{ii},'file')==0
                    [~, nrow ]=system(['cat ' obj.Params.GQI.in.bvecs{ii} ' | wc -l | awk  '' {print $1} '' '  ] );
                    nrow=str2num(nrow);
                    temp_bvecs{ii}=[ outpath 'temp.txt' ];
                    if nrow == 3 ; %then its in column form, change it...
                        exec_cmd=[ '/eris/bang/ADRC/Scripts/older/other_scripts/drigo_col2rows.sh ' obj.Params.GQI.in.bvecs{ii} ...
                            ' > ' temp_bvecs{ii}];
                        obj.RunBash(exec_cmd);
                    else
                        exec_cmd=[ 'cat ' obj.Params.GQI.in.bvecs{ii} ' >> ' temp_bvecs{ii} ];
                        obj.RunBash(exec_cmd);
                    end
                    exec_cmd=[' paste ' obj.Params.GQI.in.bvals{ii} ' ' ...
                        temp_bvecs{ii} ' | sed ''s/\t/ /g'' >' obj.Params.GQI.out.btable{ii}  ];
                    obj.RunBash(exec_cmd);
                    exec_cmd=(['rm ' temp_bvecs{ii}]);
                    obj.RunBash(exec_cmd);
                    %                     else
                    %                         fprintf(['\n B-table: ' obj.Params.GQI.out.btable{ii}  ' exists. Skipping creation...']);
                end
                
                %Attempting to create the src.fz file:
                obj.Params.GQI.out.src_fn{ii} = [ outpath obj.Params.GQI.in.prefix '.src.gz' ];
                
                if exist(obj.Params.GQI.out.src_fn{ii},'file')==0
                    fprintf('\nSource gz file reconstruction...');
                    exec_cmd=[ 'dsi_studio_run --action=src ' ...
                        ' --source=' obj.Params.GQI.in.fn{ii} ...
                        ' --b_table=' obj.Params.GQI.out.btable{ii} ...
                        ' --output=' obj.Params.GQI.out.src_fn{ii} ];
                    obj.RunBash(exec_cmd,1); %for some reason system exist with 1 :/
                    fprintf('...done\n');
                    pause(5) ; %this will add enough time for the src.gz to be completed before being read by fib.gz creation.
                    wasRun=true;
                    obj.UpdateHist(obj.Params.GQI,'proc_gqi_src', obj.Params.GQI.out.src_fn{ii},wasRun);
                else
                    [~, bb, cc ] = fileparts(obj.Params.GQI.out.src_fn{ii});
                    fprintf(['The file ' bb cc ' is complete\n']);
                end
                
                try
                    obj.Params.GQI.out.fibs_fn{ii} = ls([outpath '*.fib.gz' ] );
                catch
                    obj.Params.GQI.out.fibs_fn{ii} = '';
                end
                %Attempting to create the fib.fz file:
                if isempty(strtrim(obj.Params.GQI.out.fibs_fn{ii}))
                    fprintf('\nFib gz file reconstruction...');
                    exec_cmd=[ 'dsi_studio_run --action=rec ' ...
                        ' --source=' obj.Params.GQI.out.src_fn{ii} ...
                        ' --method=' obj.Params.GQI.in.method ...
                        ' --num_fiber=' obj.Params.GQI.in.num_fiber ...
                        ' --param0=' obj.Params.GQI.in.param0 ...
                        ' --mask=' obj.Params. GQI.in.mask{ii} ];
                    obj.RunBash(exec_cmd,1);
                    fprintf('...done\n ');
                    wasRun=true;
                    obj.UpdateHist(obj.Params.GQI,'proc_gqi_fib', strtrim(obj.Params.GQI.out.fibs_fn{ii}),wasRun);
                    
                    %Assigning the fib_fn value again (if created)
                    try
                        obj.Params.GQI.out.fibs_fn{ii} = ls([outpath '*.fib.gz' ] );
                    catch
                        obj.Params.GQI.out.fibs_fn{ii} = '';
                    end
                    
                else
                    [~, bb, cc ] = fileparts(obj.Params.GQI.out.fibs_fn{ii});
                    fprintf(['The file ' bb strtrim(cc) ' is complete \n']);
                end
                obj.Params.GQI.out.fibs_GFA{ii} = [ strtrim(obj.Params.GQI.out.fibs_fn{ii}) '.gfa.nii.gz' ];
                
                %Now exporting some values (GFA,...):
                if exist(obj.Params.GQI.out.fibs_GFA{ii},'file') == 0
                    exec_cmd=(['dsi_studio_run --action=exp ' ...
                        ' --source=' strtrim(obj.Params.GQI.out.fibs_fn{ii}) ...
                        ' --export=' obj.Params.GQI.out.export ]);
                    obj.RunBash(exec_cmd,1);
                    wasRun=true;
                    obj.UpdateHist(obj.Params.GQI,'proc_gqi_fib', strtrim(obj.Params.GQI.out.fibs_fn{ii}),wasRun);
                end
            end
        end
        
        function obj = proc_antsreg(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_ANTSREG():');
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
                    fprintf('\nCoregistering Ants to reference... ');
                    tic
                    exec_cmd=[ 'antsRegistrationSyN.sh ' ...
                        ' -d '  obj.Params.AntsReg.in.dim ...
                        ' -n '  obj.Params.AntsReg.in.threads ...
                        ' -t '  obj.Params.AntsReg.in.transform ...
                        ' -r '  obj.Params.AntsReg.in.radius  ...
                        ' -p '  obj.Params.AntsReg.in.precision ...
                        ' -f '  obj.Params.AntsReg.in.ref ...
                        ' -m '  obj.Params.AntsReg.in.fn{ii} ...
                        ' -o '  [ outpath obj.Params.AntsReg.in.prefix]  ];
                    system(exec_cmd)
                    time_taken=toc;
                    fprintf('...done.\n');
                    wasRun=true;
                    obj.UpdateHist(obj.Params.GQI,'proc_antsreg', obj.Params.AntsReg.out.fn{ii},wasRun);
                    
                else
                    [~, bb,cc ] = fileparts(obj.Params.AntsReg.out.fn{ii} );
                    fprintf(['The file ' bb cc ' is complete\n']) ;
                end
                for tocomment=1:1
                    obj.Params.AntsReg.out.FA{ii} = [ outpath obj.Params.AntsReg.in.prefix 'FA.nii.gz' ];
                    if  exist(obj.Params.AntsReg.out.FA{ii},'file')==0
                        fprintf('\n Warping dtifit metrics...');
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
                        obj.RunBash(exec_cmd);
                        %MD:
                        exec_cmd=[ 'WarpImageMultiTransform 3 ' ...
                            ' ' strrep(obj.Params.Dtifit.out.FA{ii},'FA','MD')  ...
                            ' ' strrep(obj.Params.AntsReg.out.FA{ii},'FA','MD') ...
                            ' -R '  obj.Params.AntsReg.in.ref ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped','_1Warp') ...
                            ' ' strrep(obj.Params.AntsReg.out.fn{ii},'_Warped.nii.gz','_0GenericAffine.mat') ];
                        system(exec_cmd)
                        fprintf('...done.\n');
                        obj.UpdateHist(obj.Params.GQI,'proc_antsreg_diffmetrics', strrep(obj.Params.AntsReg.out.FA{ii},'FA','RD'),wasRun);
                    end
                end
            end
        end
        
        function obj = proc_skeletonize(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_SKELETONIZE():');
            for ii=1:numel(obj.Params.Skeletonize.in.fn)
                clear cur_fn;
                if iscell(obj.Params.Skeletonize.in.fn{ii})
                    cur_fn=cell2char(obj.Params.Skeletonize.in.fn{ii});
                else
                    cur_fn=obj.Params.Skeletonize.in.fn{ii};
                end
                [a, b, c ] = fileparts(cur_fn);
                outpath=obj.getPath(a,obj.Params.Skeletonize.in.movefiles);
                clear outfile
                obj.Params.Skeletonize.out.fn{ii} = [ outpath obj.Params.Skeletonize.in.prefix '.nii.gz' ];
                if exist(obj.Params.Skeletonize.out.fn{ii},'file')==0
                    fprintf('\nSkeletonizing to reference... ');
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
                    fprintf('...done.\n');
                    wasRun=true;
                    obj.UpdateHist(obj.Params.Skeletonize,'proc_skeletonize', obj.Params.Skeletonize.out.fn{ii},wasRun);
                else
                    [ ~ , bb, cc ] = fileparts(obj.Params.AntsReg.out.fn{ii});
                    fprintf(['The file ' bb cc ' is comeplete.\n']) ;
                end
                
                
                %NOW OTHER METRICS:
                obj.Params.Skeletonize.in.FA{ii}  = obj.Params.AntsReg.out.FA{ii} ;
                obj.Params.Skeletonize.out.FA{ii} = [ outpath obj.Params.Skeletonize.in.prefix '_FA.nii.gz' ];
                
                for tocomment=1:1;
                    if exist(obj.Params.Skeletonize.out.FA{ii},'file')==0
                        fprintf('\nSkeletonizing  dtimetrics... ');
                        tic
                        fprintf('\n in FA...');
                        %For FA we move the projected skeletonize into FA
                        %(as we use this metric for the initial guess)
                        exec_cmd =[ 'cp ' obj.Params.Skeletonize.out.fn{ii} ' ' obj.Params.Skeletonize.out.FA{ii} ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done\n');
                        
                        %For RD/AxD//MD a -a (alternative metrics) should
                        %be added...
                        %RD:
                        fprintf('\n in RD...');
                        exec_cmd=[ 'tbss_skeleton ' ...
                            ' -i '  obj.Params.Skeletonize.in.meanFA ...
                            ' -p '  obj.Params.Skeletonize.in.thr ...
                            ' '  obj.Params.Skeletonize.in.skel_dst ...
                            ' '  obj.Params.Skeletonize.in.ref_region  ...
                            ' '  obj.Params.Skeletonize.in.FA{ii} ...
                            ' '  strrep(obj.Params.Skeletonize.out.FA{ii},'FA','RD') ...
                            ' -a ' strrep(obj.Params.Skeletonize.in.FA{ii},'FA','RD')  ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done\n');
                        %AxD:
                        fprintf('\n in AxD...');
                        exec_cmd=[ 'tbss_skeleton ' ...
                            ' -i '  obj.Params.Skeletonize.in.meanFA ...
                            ' -p '  obj.Params.Skeletonize.in.thr ...
                            ' '  obj.Params.Skeletonize.in.skel_dst ...
                            ' '  obj.Params.Skeletonize.in.ref_region  ...
                            ' '  obj.Params.Skeletonize.in.FA{ii} ...
                            ' '  strrep(obj.Params.Skeletonize.out.FA{ii},'FA','AxD') ...
                            ' -a ' strrep(obj.Params.Skeletonize.in.FA{ii},'FA','AxD')  ];
                        
                        
                        obj.RunBash(exec_cmd); fprintf('...done\n');
                        %MD:
                        fprintf('\n in MD...');
                        exec_cmd=[ 'tbss_skeleton ' ...
                            ' -i '  obj.Params.Skeletonize.in.meanFA ...
                            ' -p '  obj.Params.Skeletonize.in.thr ...
                            ' '  obj.Params.Skeletonize.in.skel_dst ...
                            ' '  obj.Params.Skeletonize.in.ref_region  ...
                            ' '  obj.Params.Skeletonize.in.FA{ii} ...
                            ' '  strrep(obj.Params.Skeletonize.out.FA{ii},'FA','MD') ...
                            ' -a ' strrep(obj.Params.Skeletonize.in.FA{ii},'FA','MD')  ];
                        
                        obj.RunBash(exec_cmd); fprintf('...done\n');
                        toc
                        fprintf('...done.\n');
                    else
                        [~, bb, cc ] = fileparts(obj.Params.AntsReg.out.fn{ii});
                        fprintf(['The file ' bb cc ' is complete.\n']) ;
                        
                    end
                end
                obj.Params.Skeletonize.out.diffmetrics={ 'FA' 'RD' 'AxD' 'MD' } ;
            end
        end
        
        function obj = proc_getskeltois(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_GETSKELTOIS():');
            for kk=1:numel(obj.Params.Skeletonize.out.FA)
                for jj=1:numel( obj.Params.Skeletonize.out.diffmetrics)
                    for ii=1:numel(obj.Params.Skel_TOI.in.masks)
                        cur_name = [ obj.Params.Skeletonize.out.diffmetrics{jj} '_' obj.Params.Skel_TOI.in.masks{ii} obj.Params.Skel_TOI.in.suffix ] ;
                        if ~isfield(obj.Params.Skel_TOI.out,cur_name)
                            obj.Params.Skel_TOI.out.(cur_name) = '';
                        end
                        if isempty(obj.Params.Skel_TOI.out.(cur_name)) || size(obj.Params.Skel_TOI.out.(cur_name),2) ~= 10 %if not 10 characters, thena mistake occured!
                            cur_field=[ obj.Params.Skeletonize.out.diffmetrics{jj} ...
                                '_' obj.Params.Skel_TOI.in.masks{ii}  obj.Params.Skel_TOI.in.suffix ];
                            in_file=strrep(obj.Params.Skeletonize.out.FA{kk},'_FA.nii',[ '_' obj.Params.Skeletonize.out.diffmetrics{jj} '.nii' ] );
                            mask_file=[ obj.Params.Skel_TOI.in.location obj.Params.Skel_TOI.in.masks{ii}   '.nii.gz' ] ;
                            
                            exec_cmd=['fslstats ' in_file ' -k ' mask_file ' -M '  ];
                            fprintf([ ' now in ' cur_name '\n'] );
                            [~ , obj.Params.Skel_TOI.out.(cur_name) ] =  system(exec_cmd);
                            last_cur_name=cur_name;
                            wasRun=true;
                        end
                        clear cur_field  in_file out_file mask_file cur_name ;
                    end
                end
                wasRun=true;
                obj.UpdateHist(obj.Params.Skel_TOI,'proc_getskeltois', '',wasRun);
                clear last_cur_name;
                fprintf('...done\n');
            end
            fprintf('proc_getskeltois() is complete.\n')
        end
        
        function obj = proc_getFreeSurfer(obj)
            wasRun=false;
            %Double check what is your default shell (to source
            %FreeSurfer):;;
            %if strcmp(obj.Params.FreeSurfer.shell,'bash')
            fprintf('\n%s\n', 'PERFORMING PROC_GETFREESURFER():');
            %display([ 'The whoami output is: ' obj.Params.FreeSurfer.shell ])
            
            if exist(obj.Params.FreeSurfer.out.aparcaseg, 'file') == 0
                if strcmp(obj.Params.FreeSurfer.shell,'rdp20') %due to launchpad errors, I decided to use this 'whoami' instead of shell. NEED TO FIX IT!
                    export_shell=[ 'export FREESURFER_HOME=' obj.Params.FreeSurfer.init_location ' ; '...
                        ' source $FREESURFER_HOME/SetUpFreeSurfer.sh ;' ...
                        ' export SUBJECTS_DIR=' obj.Params.FreeSurfer.dir ' ; '];
                else
                    export_shell=[ ' setenv FREESURFER_HOME ' obj.Params.FreeSurfer.init_location ' ; ' ...
                        ' source $FREESURFER_HOME/SetUpFreeSurfer.csh ; ' ...
                        ' setenv SUBJECTS_DIR ' obj.Params.FreeSurfer.dir ' ; '];
                end
                
                
                %Attempting to create B0means:
                if obj.Params.FreeSurfer.in.T2exist==true
                    %use T2 for recon-all
                    exec_cmd=[ export_shell ...
                        ' recon-all -all -subjid ' obj.sessionname ...
                        ' -deface -i ' strtrim(obj.Params.FreeSurfer.in.T1) ...
                        ' -T2 ' strtrim(obj.Params.FreeSurfer.in.T2) ...
                        ' -hippocampal-subfields-T1T2 ' strtrim(obj.Params.FreeSurfer.in.T2) ' T2 ' ];
                else
                    %no T2 so use only T1 for recon-all
                    exec_cmd=[ export_shell ...
                        ' recon-all -all  -subjid ' obj.sessionname ...
                        ' -deface -i ' strtrim(obj.Params.FreeSurfer.in.T1) ...
                        ' -hippocampal-subfields-T1' ];
                end
                disp('Running FreeSurfer... (this will take ~24 hours)')
                tic;
                obj.RunBash(exec_cmd,44); % '44' codes for seeing the output!
                obj.Params.FreeSurfer.out.timelapsed_mins=toc/60;
                disp('Done with FreeSurfer');
                wasRun=true;
                obj.UpdateHist(obj.Params.FreeSurfer,'proc_getFS', obj.Params.FreeSurfer.out.aparcaseg,wasRun);
            else
                [~,bb,cc ] = fileparts(obj.Params.FreeSurfer.out.aparcaseg) ;
                fprintf(['The aparc aseg file ' bb cc ' exists. \n' ]);
            end
        end
        
        function obj = proc_FS2dwi(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_FS2DWI():');
            
            %INIT SPECS:
            for tohide=1:1
                %Create folder to put ouput:
                [a, ~, ~] = fileparts(obj.Params.FS2dwi.in.b0{1});
                outpath=obj.getPath(a,obj.Params.FS2dwi.in.movefiles );
                %Init outputs:
                obj.Params.FS2dwi.out.xfm_dwi2FS = [ outpath 'xfm_dwi2FS.lta' ] ;
                obj.Params.FS2dwi.out.fn_aparc = [ outpath  'dwi_aparc+aseg.nii.gz' ] ;
                obj.Params.FS2dwi.out.fn_aparc2009 = [ outpath 'dwi_aparc.a2009+aseg.nii.gz' ] ;
                obj.Params.FS2dwi.out.hippofield_left = [ outpath  'dwi_hippofields_lh.nii.gz' ] ;
                obj.Params.FS2dwi.out.hippofield_right = [ outpath  'dwi_hippofields_rh.nii.gz' ] ;
                
                %Making directory for all the other output
                system (['mkdir -p ' outpath filesep 'aparc_aseg' filesep]);
                system (['mkdir -p ' outpath filesep 'aparc2009_aseg' filesep]);
                system (['mkdir -p ' outpath filesep 'hippos' filesep]);
                
                %Sourcing FS and SUBJECTS_DIR:
                if strcmp(obj.Params.FreeSurfer.shell,'rdp20') %due to launchpad errors, I decided to use this 'whoami' instead of shell. NEED TO FIX IT!
                    export_shell=[ 'export FREESURFER_HOME=' obj.Params.FreeSurfer.init_location ' ; '...
                        ' source $FREESURFER_HOME/SetUpFreeSurfer.sh ;' ...
                        ' export SUBJECTS_DIR=' obj.Params.FreeSurfer.dir ' ; '];
                else
                    export_shell=[ ' setenv FREESURFER_HOME ' obj.Params.FreeSurfer.init_location ' ; ' ...
                        ' source $FREESURFER_HOME/SetUpFreeSurfer.csh ; ' ...
                        ' setenv SUBJECTS_DIR ' obj.Params.FreeSurfer.dir ' ; '];
                end
            end
            %APARCASEG EXTRACTION:
            for tohide=1:1
                if exist(obj.Params.FS2dwi.in.aparcaseg, 'file') == 2
                    %BBreg dwi (b0) to FS:
                    if exist(obj.Params.FS2dwi.out.xfm_dwi2FS,'file') == 0
                        %bbreg b0 to FS_T1:
                        exec_cmd=[ export_shell ...
                            ' bbregister --s ' obj.sessionname ...
                            ' --mov ' obj.Params.FS2dwi.in.b0{1} ...
                            ' --reg ' obj.Params.FS2dwi.out.xfm_dwi2FS ' --dti --init-fsl '];
                        disp('proc_FS2dwi: Running bbreg dwi2FS_T1...')
                        obj.RunBash(exec_cmd,44); % '44' codes for seeing the output!
                        disp('..done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.FS2dwi,'proc_FS2dwi', obj.Params.FS2dwi.out.xfm_dwi2FS, wasRun);
                    end
                    
                    %Aparc+aseg to dwi:
                    if exist(obj.Params.FS2dwi.out.fn_aparc,'file') == 0
                        %bbreg b0 to FS_T1:
                        exec_cmd=[ export_shell ' mri_vol2vol --mov ' obj.Params.FS2dwi.in.b0{1} ...
                            ' --targ ' obj.Params.FS2dwi.in.aparcaseg ...
                            ' --o ' obj.Params.FS2dwi.out.fn_aparc ...
                            ' --inv --nearest --reg ' obj.Params.FS2dwi.out.xfm_dwi2FS  ];
                        disp('proc_FS2dwi: Running bbreg in aparc+aseg.mgz...')
                        obj.RunBash(exec_cmd,44); % '44' codes for seeing the output!
                        disp('..done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.FS2dwi,'proc_FS2dwi', obj.Params.FS2dwi.out.fn_aparc , wasRun);
                    end
                    
                    %Extracting all ROIs for aparc:
                    [num_aparc name_aparc ] =textread(obj.Params.FS2dwi.in.tmpfile_aparcaseg,'%s %s');
                    for ff=1:numel(num_aparc)
                        tmp_curname{ff} = [ outpath  'aparc_aseg' filesep 'dwi_' name_aparc{ff} '.nii.gz'];
                        if exist(strtrim(tmp_curname{ff}), 'file') == 0
                            fprintf(['\nDisplaying now: ' tmp_curname{ff} '...' ] )
                            exec_cmd = [ 'fslmaths  ' obj.Params.FS2dwi.out.fn_aparc ...
                                ' -uthr ' num_aparc{ff} ' -thr ' num_aparc{ff} ...
                                ' -div '  num_aparc{ff} ' ' tmp_curname{ff} ] ;
                            fprintf('done \n')
                            obj.RunBash(exec_cmd);
                        end
                    end
                    fprintf('aparc+aseg extraction complete\n')
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    error([ 'in proc_FS2dwi(obj): FS aparc_aseg (supposely) located in: ' obj.Params.FS2dwi.in.aparcaseg ' does not exist'])
                end
                fprintf('aparc+aseg extraction complete. \n')
            end
            %APARCASEGa2009 EXTRACTION:
            for tohide=1:1
                if exist(obj.Params.FS2dwi.in.aparcaseg2009,'file')==2
                    %Aparc.a2009+aseg to dwi:
                    if exist(obj.Params.FS2dwi.out.fn_aparc2009,'file') == 0
                        %bbreg b0 to FS_T1:
                        exec_cmd=[ export_shell ' mri_vol2vol --mov ' obj.Params.FS2dwi.in.b0{1} ...
                            ' --targ ' obj.Params.FS2dwi.in.aparcaseg2009 ...
                            ' --o '  obj.Params.FS2dwi.out.fn_aparc2009 ...
                            ' --inv --nearest --reg ' obj.Params.FS2dwi.out.xfm_dwi2FS  ];
                        disp('proc_FS2dwi: Running bbreg in aparc2009+aseg.mgz...')
                        obj.RunBash(exec_cmd,44); % '44' codes for seeing the output!
                        disp('..done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.FS2dwi,'proc_FS2dwi', obj.Params.FS2dwi.out.fn_aparc2009 , wasRun);
                    end
                    %Extracting all ROIs for aparc2009:
                    [num_aparc2009 name_aparc2009 ] =textread(obj.Params.FS2dwi.in.tmpfile_aparcaseg2009,'%s %s');
                    clear tmp_curname;
                    for ff=1:numel(num_aparc2009)
                        tmp_curname{ff} = [ outpath  'aparc2009_aseg' filesep  'dwi_' name_aparc2009{ff}  '.nii.gz' ];
                        if exist(strtrim(tmp_curname{ff}), 'file') == 0
                            fprintf(['\nDisplaying now: ' tmp_curname{ff} '...' ] )
                            exec_cmd = [ 'fslmaths  ' obj.Params.FS2dwi.out.fn_aparc ...
                                ' -uthr ' num_aparc2009{ff} ' -thr ' num_aparc2009{ff} ...
                                ' -div '  num_aparc2009{ff} ' ' tmp_curname{ff} ] ;
                            fprintf('done \n')
                            obj.RunBash(exec_cmd);
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    error([ 'in proc_FS2dwi(obj): FS aparc_aseg (supposely) located in: ' obj.Params.FS2dwi.out.fn_aparc2009 ' does not exist'])
                end
                fprintf('aparc+aseg2009 extraction complete. \n')
            end
            %HIPPOFIELD_LEFT EXTRACTION:
            for tohide=1:1
                if exist( obj.Params.FS2dwi.in.hippofield_left ,'file') == 2
                    if exist(obj.Params.FS2dwi.out.hippofield_left,'file') == 0
                        %bbreg b0 to FS_T1:
                        exec_cmd=[ ' mri_vol2vol --mov ' obj.Params.FS2dwi.in.b0{1} ...
                            ' --targ ' obj.Params.FS2dwi.in.hippofield_left ...
                            ' --o ' obj.Params.FS2dwi.out.hippofield_left ...
                            ' --inv --nearest --reg ' obj.Params.FS2dwi.out.xfm_dwi2FS  ];
                        disp('proc_FS2dwi: Running bbreg in aparc+aseg.mgz...')
                        obj.RunBash(exec_cmd); % '44' codes for seeing the output!
                        disp('..done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.FS2dwi,'proc_FS2dwi', obj.Params.FS2dwi.out.hippofield_left , wasRun);
                    end
                    [num_hippo_lh name_hippo_lh ] =textread(obj.Params.FS2dwi.in.tmpfile_hippo_bil,'%s %s');
                    clear tmp_curname;
                    for ff=1:numel(num_hippo_lh)
                        tmp_curname{ff} = [ outpath  'hippos' filesep 'dwi_lh_' name_hippo_lh{ff} '.nii.gz' ];
                        if exist(strtrim(tmp_curname{ff}), 'file') == 0
                            fprintf(['\nDisplaying now: ' tmp_curname{ff} '...' ] )
                            exec_cmd = [ 'fslmaths  ' obj.Params.FS2dwi.out.hippofield_left ...
                                ' -uthr ' num_hippo_lh{ff} ' -thr ' num_hippo_lh{ff} ...
                                ' -div '  num_hippo_lh{ff} ' ' tmp_curname{ff} ] ;
                            fprintf('done \n')
                            obj.RunBash(exec_cmd);
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    error([ 'in proc_FS2dwi(obj): FS hippofield_left (supposely) located in: ' obj.Params.FS2dwi.in.hippofield_left ' does not exist'])
                end
                fprintf('hippofield_lh extraction complete\n')
            end
            
            %HIPPOFIELD_RIGHT EXTRACTION:
            for tohide=1:1
                if exist( obj.Params.FS2dwi.in.hippofield_right ,'file') == 2
                    if exist(obj.Params.FS2dwi.out.hippofield_right,'file') == 0
                        %bbreg b0 to FS_T1:
                        exec_cmd=[ ' mri_vol2vol --mov ' obj.Params.FS2dwi.in.b0{1} ...
                            ' --targ ' obj.Params.FS2dwi.in.hippofield_right ...
                            ' --o ' obj.Params.FS2dwi.out.hippofield_right ...
                            ' --inv --nearest --reg ' obj.Params.FS2dwi.out.xfm_dwi2FS  ];
                        disp('proc_FS2dwi: Running bbreg in aparc+aseg.mgz...')
                        obj.RunBash(exec_cmd); % '44' codes for seeing the output!
                        disp('..done');
                        wasRun=true;
                        obj.UpdateHist(obj.Params.FS2dwi,'proc_FS2dwi', obj.Params.FS2dwi.out.hippofield_right , wasRun);
                    end
                    
                    %Extracting all ROIs for aparc2009:
                    [num_hippo_rh name_hippo_rh ] =textread(obj.Params.FS2dwi.in.tmpfile_hippo_bil,'%s %s');
                    clear tmp_curname;
                    for ff=1:numel(num_hippo_rh)
                        tmp_curname{ff} = [outpath  'hippos' filesep 'dwi_rh_' name_hippo_rh{ff} '.nii.gz' ];
                        if exist(strtrim(tmp_curname{ff}), 'file') == 0
                            fprintf(['\nDisplaying now: ' tmp_curname{ff} '...' ] )
                            exec_cmd = [ 'fslmaths  ' obj.Params.FS2dwi.out.hippofield_right ...
                                ' -uthr ' num_hippo_rh{ff} ' -thr ' num_hippo_rh{ff} ...
                                ' -div '  num_hippo_rh{ff} ' ' tmp_curname{ff} ] ;
                            fprintf('done \n')
                            obj.RunBash(exec_cmd);
                        end
                    end
                else
                    error([ 'in proc_FS2dwi(obj): FS hippofield_right (supposely) located in: ' obj.Params.FS2dwi.in.hippofield_right ' does not exist'])
                end
                fprintf('hippofield_rh extraction complete\n')
            end
        end
        
        function obj = proc_FROIS2dwi(obj)
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_FROIS2DWI():');
            [a, ~, ~ ] = fileparts(obj.Params.FROIS2dwi.in.fn);
            outpath=obj.getPath(a,obj.Params.FROIS2dwi.in.movefiles);
            
            %make a list of all *.gz that exist in that direcotry
            [ok_check, tmp_frois_list ]  = system(['ls -1 '  ...
                obj.Params.FROIS2dwi.in.FROIS_dir '*.nii.gz']);
            if ok_check ~= 0
                error(['In proc_FROIS2dwi(): cannot find any *.nii.gz images in:'...
                    obj.Params.FROIS2dwi.in.FROIS_dir 'Double check']);
            end
            
            tmp_cellarray=textscan(tmp_frois_list,'%s');
            obj.Params.FROIS2dwi.in.FROIS_list=tmp_cellarray{1};
            
            
            %Check if the MNI_T1 exists....
            if exist(obj.Params.FROIS2dwi.in.MNI_T1,'file') == 0
                error(['In proc_FROIS2dwi(): Cannnot find the MNI_T1 in:' obj.Params.FROIS2dwi.in.MNI_T1 ] );
            end
            
            %Assigend outpath of coreg using ants
            obj.Params.FROIS2dwi.out.MNI_2_dwi = [ outpath obj.Params.FROIS2dwi.in.prefix 'Warped.nii.gz' ];
            if exist(obj.Params.FROIS2dwi.out.MNI_2_dwi ,'file')==0
                fprintf('\n In proc_FROIS(): Coregistering Ants to reference... ');
                aa=1;
                %obj.proc_as_coreg(obj.Params.FROIS2dwi.in.FROIS_list
                fprintf('...done.\n');
            end
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%% END Data Pre-Processing Methods %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% BEGIN Data Post-Processing Methods %%%%%%%%%%%%%
        function obj = proc_tracula(obj)
            % try
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_TRACULA():');
            %Creating root directory:
            for tohide=1:1
                [a b c ] = fileparts(obj.Params.Tracula.in.fn);
                outpath=obj.getPath(a,obj.Params.Tracula.in.movefiles);
                obj.Params.Tracula.out.dcmirc = [outpath 'dcmrirc.' obj.Params.Tracula.in.prefix] ;
                
                if strcmp(obj.projectID,'HAB')
                    replaced_outpath = outpath ;
                    outpath = [ '/eris/bang/HAB_Project1/TRACULA' filesep obj.sessionname filesep ];
                    system(['mkdir -p ' outpath ]);
                    system(['ln -s ' outpath ' ' replaced_outpath filesep obj.sessionname ]);
                end
            end
            %Create the necessary dcmirc file:
            for tohide=1:1
                if exist(obj.Params.Tracula.out.dcmirc,'file') == 0
                    exec_cmd = [ 'cat ' obj.Params.Tracula.in.dcmrirc ' | sed s%''<SUBJECTID>''%' ...
                        obj.sessionname '%g | sed s%''<DTROOT>''%' outpath ...
                        '%g | sed s%''<FSSUBJECTSDIR>''%' obj.Params.Tracula.in.FSDIR ...
                        '%g | sed s%''<BVECFILE>''%' obj.Params.Tracula.in.bvec ...
                        '%g | sed s%''<BVALFILE>''%' obj.Params.Tracula.in.bval ...
                        '%g | sed s%''<NB0>''%' num2str(obj.Params.Tracula.in.nb0) ...
                        '%g > ' obj.Params.Tracula.out.dcmirc ];
                    obj.RunBash(exec_cmd);
                else
                    [~, bb, cc ] = fileparts(obj.Params.Tracula.out.dcmirc);
                    fprintf(['The file ' bb cc ' exists.\n']);
                end
            end
            %Run the three necessary steps for TRACULA (and bedpostx)
            for tohide=1:1
                %Set output directory:
                obj.Params.Tracula.out.dir = outpath ;
                %
                %Check for final file for step 1:;;;;
                obj.Params.Tracula.out.prep_check = [ obj.Params.Tracula.out.dir  obj.sessionname ...
                    filesep 'dlabel' filesep 'diff' filesep 'lh.cst_AS_avg33_mni_bbr_cpts_6.nii.gz' ];
                
                obj.Params.Tracula.out.isrunning = [ obj.Params.Tracula.out.dir  obj.sessionname ...
                    filesep  'scripts' filesep 'IsRunning.trac' ];
                
                if exist(obj.Params.Tracula.out.prep_check,'file') == 0
                    %Remove IsRunning.trac if exists
                    if exist( obj.Params.Tracula.out.isrunning, 'file') ~= 0
                        system(['rm '  obj.Params.Tracula.out.isrunning ]);
                    end
                    exec_cmd = ['trac-all -prep -c ' obj.Params.Tracula.out.dcmirc ' -i ' obj.Params.Tracula.in.fn ];
                    obj.RunBash(exec_cmd,44);
                else
                    [~, bb, cc ] = fileparts(obj.Params.Tracula.out.prep_check);
                    fprintf(['trac-all -prep filecheck ' bb cc ' exists.\n']);
                end
                %Step 2: Trac-all -prep
                obj.Params.Tracula.out.bedp_check = [ obj.Params.Tracula.out.dir  obj.sessionname ...
                    filesep 'dmri.bedpostX' filesep 'mean_fsumsamples.nii.gz' ];
                if exist(obj.Params.Tracula.out.bedp_check, 'file') == 0
                    %Remove IsRunning.trac if exists
                    if exist( obj.Params.Tracula.out.isrunning, 'file') ~= 0
                        system(['rm '  obj.Params.Tracula.out.isrunning ]);
                    end
                    exec_cmd = ['trac-all -bedp -c ' obj.Params.Tracula.out.dcmirc ' -i ' obj.Params.Tracula.in.fn ];
                    obj.RunBash(exec_cmd,44);
                else
                    [~, bb, cc ] = fileparts(obj.Params.Tracula.out.bedp_check);
                    fprintf(['trac-all -bedp file ' bb cc ' exists.\n']);
                end
                %Step 3: Trac-all -path
                obj.Params.Tracula.out.path_check = [ obj.Params.Tracula.out.dir  obj.sessionname ...
                    filesep 'dpath' filesep 'merged_avg33_mni_bbr.mgz' ];
                if exist(obj.Params.Tracula.out.path_check,'file') == 0
                    %Remove IsRunning.trac if exists
                    if exist( obj.Params.Tracula.out.isrunning, 'file') ~= 0
                        system(['rm '  obj.Params.Tracula.out.isrunning ]);
                    end
                    exec_cmd = ['trac-all -path -c ' obj.Params.Tracula.out.dcmirc ' -i ' obj.Params.Tracula.in.fn ];
                    obj.RunBash(exec_cmd,44);
                else
                    [~, bb, cc ] = fileparts(obj.Params.Tracula.out.path_check);
                    fprintf(['trac-all -path ' bb cc ' exists.\n']);
                end
                
            end
        end
        
        function obj = proc_AFQ(obj)
            % try
            wasRun=false;
            fprintf('\n%s\n', 'PERFORMING PROC_AFQ():');
            %Creating root directory:
            for tohide=1:1
                [a b c ] = fileparts(obj.Params.AFQ.in.dwi);
                outpath=obj.getPath(a,obj.Params.AFQ.in.movefiles);
            end
            %Run the three necessary steps for AFQ
            for tohide=1:1
                %Set output directory:
                obj.Params.AFQ.out.dir = outpath ;
                obj.Params.AFQ.out.dwi = 'AFQ_dn.mat';
                %
                [obj.Params.AFQ.out.dwi, obj.Params.AFQ.out.dir] = dtiInit(obj.Params.AFQ.in.dwi, obj.Params.AFQ.in.T1, []);
                
                
            end
        end
        
        function obj = trkland_fx(obj)
            wasRun = false;
            fprintf('\n%s\n', 'PERFORMING TRKLAND FORNIX: TRKLAND_FX():');
            %Initialize which image you'll use for FLIRT (either FA or b0)
            %This is changeable!
            
            obj.Trkland.fx.in.ref = obj.Trkland.fx.in.b0;
            %obj.Trkland.fx.in.ref = obj.Trkland.fx.in.FA;
            
            %Creating root directory:
            for tohide=1:1
                [a, b, c ] = fileparts(obj.Trkland.fx.in.b0);
                outpath=obj.getPath(a,obj.Trkland.fx.in.movefiles);
                obj.Trkland.fx.out.QCfile_lh = [outpath 'QC_fx_lh.flag'] ;
                obj.Trkland.fx.out.QCfile_rh = [outpath 'QC_fx_rh.flag'] ;
                obj.Trkland.fx.out.QCfile_bil = [outpath 'QC_fx_bil.flag'] ;
                
            end
            
            
            %Remove a field that is not being used anymore:
            if isfield(obj.Trkland.fx.out,'QC')
                rmfield(obj.Trkland.fx.out,'QC');
            end
            %MATFILE TRANSFORMATION SECTION:
            for tohide=1:1
                %Create fx directory (if doesn't exist)
                exec_cmd = [ 'mkdir -p ' obj.Trkland.root ];
                obj.RunBash(exec_cmd);
                %Create the matfile of the fx tranformation
                if exist(obj.Trkland.fx.in.tmp2b0_matfile,'file') == 0
                    %USING FLIRT, HERE WE COULD ALSO TRY SPM_COREG (not sure if
                    %I'll get/how to get the *.mat though...). Flirt does a decent job anyhow
                    fprintf(['\n Coregistering TMP-B0 to:' obj.Trkland.fx.in.ref ] );
                    exec_cmd = ['flirt -in ' obj.Trkland.fx.tmp.b0  ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -omat ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -out ' obj.Trkland.fx.in.fn_tmp2b0 ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
            %%%%%%%%%%%%%%%%%%%%%
            %BILATERAL SECTION:
            for tohide=1:1
                %Apply the matfile to the roi:
                if exist(obj.Trkland.fx.in.roi_bil, 'file') == 0
                    fprintf('\n Coregistering bil_roi...')
                    exec_cmd = ['flirt -in '  obj.Trkland.fx.tmp.roi_bil  ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roi_bil ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Apply the matfile to the dilated solid fx:
                if exist(obj.Trkland.fx.in.roa_bil_solid, 'file') == 0
                    fprintf('\n Coregistering bil_solid...')
                    exec_cmd = ['flirt -in ' obj.Trkland.fx.tmp.roa_solid_bil ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roa_bil_solid ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Erode the solid ROA:
                if exist(obj.Trkland.fx.in.roa_bil_ero , 'file') == 0
                    fprintf('\n Eroding TMP_B0_fx_solid...')
                    exec_cmd = ['fslmaths ' obj.Trkland.fx.in.roa_bil_solid ...
                        ' -ero ' obj.Trkland.fx.in.roa_bil_ero  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Creating the hollow ROA:
                if exist(obj.Trkland.fx.in.roa_bil_hollow, 'file') == 0
                    fprintf('\n Hollowing the ROA...')
                    exec_cmd = ['fslmaths ' obj.Trkland.fx.in.roa_bil_solid ...
                        ' -sub ' obj.Trkland.fx.in.roa_bil_ero ' ' obj.Trkland.fx.in.roa_bil_hollow ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
            %%%%%%%%%%%%%%%%%%
            %LEFT HEMISPHERE SECTION:
            for tohide=1:1
                %Apply the matfile to the roi:
                if exist(obj.Trkland.fx.in.roi_lh, 'file') == 0
                    fprintf('\n Coregistering lh_roi...')
                    exec_cmd = ['flirt -in '  obj.Trkland.fx.tmp.roi_lh  ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roi_lh ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Apply the matfile to the dilated solid fx:
                if exist(obj.Trkland.fx.in.roa_lh_solid, 'file') == 0
                    fprintf('\n Coregistering lh_solid...')
                    exec_cmd = ['flirt -in ' obj.Trkland.fx.tmp.roa_solid_lh ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roa_lh_solid ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Erode the solid ROA:
                if exist(obj.Trkland.fx.in.roa_lh_ero , 'file') == 0
                    fprintf('\n Eroding TMP_B0_fx_solid...')
                    exec_cmd = ['fslmaths ' obj.Trkland.fx.in.roa_lh_solid ...
                        ' -ero ' obj.Trkland.fx.in.roa_lh_ero  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Creating the hollow ROA:
                if exist(obj.Trkland.fx.in.roa_lh_hollow, 'file') == 0
                    fprintf('\n Hollowing the ROA...')
                    exec_cmd = ['fslmaths ' obj.Trkland.fx.in.roa_lh_solid ...
                        ' -sub ' obj.Trkland.fx.in.roa_lh_ero ' ' obj.Trkland.fx.in.roa_lh_hollow ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
            %%%%%%%%%%%%%%%%%%
            %RIGHT HEMISPHERE SECTION:
            for tohide=1:1
                %Apply the matfile to the roi:
                if exist(obj.Trkland.fx.in.roi_rh, 'file') == 0
                    fprintf('\n Coregistering rh_roi...')
                    exec_cmd = ['flirt -in '  obj.Trkland.fx.tmp.roi_rh  ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roi_rh ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Apply the matfile to the dilated solid fx:
                if exist(obj.Trkland.fx.in.roa_rh_solid, 'file') == 0
                    fprintf('\n Coregistering rh_solid...')
                    exec_cmd = ['flirt -in ' obj.Trkland.fx.tmp.roa_solid_rh ...
                        ' -ref ' obj.Trkland.fx.in.ref ' -applyxfm -init ' obj.Trkland.fx.in.tmp2b0_matfile ...
                        ' -interp  nearestneighbour -out ' obj.Trkland.fx.in.roa_rh_solid ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Erode the solid ROA:
                if exist(obj.Trkland.fx.in.roa_rh_ero , 'file') == 0
                    fprintf('\n Eroding TMP_B0_fx_solid...')
                    exec_cmd = ['fslmaths '  obj.Trkland.fx.in.roa_rh_solid ...
                        ' -ero ' obj.Trkland.fx.in.roa_rh_ero  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %Creating the hollow ROA:
                if exist(obj.Trkland.fx.in.roa_rh_hollow, 'file') == 0
                    fprintf('\n Hollowing the ROA...')
                    exec_cmd = ['fslmaths ' obj.Trkland.fx.in.roa_rh_solid ...
                        ' -sub ' obj.Trkland.fx.in.roa_rh_ero ' ' obj.Trkland.fx.in.roa_rh_hollow ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
            
            %%%%%%%%%%%%%%%%%%
            %INIT TRIMMED OUTPUTS
            for tohide=1:1
                
                %INIT OUTPUTS
                obj.Trkland.fx.out.clean_trks_lh = [ obj.Trkland.root  'trkk_fx_trimmedclean_lh.trk.gz'];
                obj.Trkland.fx.out.clean_trks_rh = [ obj.Trkland.root  'trkk_fx_trimmedclean_rh.trk.gz'];
                obj.Trkland.fx.out.clean_trkstrimmed_lh = [  obj.Trkland.root  'trkk_fx_trimmed_lh.trk.gz' ];
                obj.Trkland.fx.out.clean_trkstrimmed_rh = [  obj.Trkland.root  'trkk_fx_trimmed_rh.trk.gz'];
                obj.Trkland.fx.out.clineFA_lh_highFA = [ obj.Trkland.root  'cline_fx_highFA_lh.trk.gz'];
                obj.Trkland.fx.out.clineFA_lh_HDorff = [ obj.Trkland.root  'cline_fx_HDorff_lh.trk.gz'];
                obj.Trkland.fx.out.clineFA_rh_highFA = [ obj.Trkland.root  'cline_fx_highFA_rh.trk.gz'];
                obj.Trkland.fx.out.clineFA_rh_HDorff = [ obj.Trkland.root  'cline_fx_HDorff_rh.trk.gz'];
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STARTING THE ACTUAL TRACKING NOW:
            for tohide=1:1
                %Init outputs:
                obj.Trkland.fx.out.trks_lh = [ obj.Trkland.root  'trkk_fx_lh.trk.gz'];
                obj.Trkland.fx.out.trks_rh = [ obj.Trkland.root  'trkk_fx_rh.trk.gz'];
                
                if exist(obj.Trkland.fx.out.QCfile_bil,'file') == 0
                    if exist(obj.Trkland.fx.out.QCfile_lh) == 0
                        %Left side trking:
                        if exist(obj.Trkland.fx.out.trks_lh,'file') == 0
                            if strcmp(obj.projectID,'ADRC')
                                exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                    ' --seed_count=10000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                    ' --seed=' obj.Trkland.fx.in.roi_lh ' --roa=' obj.Trkland.fx.in.roa_lh_hollow ...
                                    ' --threshold_index=nqa  --fa_threshold=0.04 --fiber_count=500' ...
                                    ' --step_size=1 --turning_angle=40 --min_length=80 --max_length=250 ' ...
                                    ' --output=' obj.Trkland.fx.out.trks_lh ];
                            else
                                exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                    ' --seed_count=10000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                    ' --seed=' obj.Trkland.fx.in.roi_lh ' --roa=' obj.Trkland.fx.in.roa_lh_hollow ...
                                    ' --step_size=1 --turning_angle=40 --min_length=80 --max_length=250 ' ...
                                    ' --output=' obj.Trkland.fx.out.trks_lh ];
                            end
                            for dd=1:4 %trying 4 times to get a trk. If not, quit!
                                if exist(obj.Trkland.fx.out.trks_lh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            wasRun=true;
                            obj.UpdateHist(obj.Trkland.fx,'trkland_fx', obj.Trkland.fx.out.trks_lh,wasRun);
                        else
                            [~, bb, cc ] = fileparts(obj.Trkland.fx.in.roa_lh_hollow);
                            fprintf(['The file ' bb cc ' exists. \n']);
                        end
                    else
                        display('QC_flag_lh found in trklnad_fx. Skipping and removing data points...')
                        %RefreshFields(obj,'fx','lh')
                    end
                    
                    if exist(obj.Trkland.fx.out.QCfile_rh) == 0
                        %Right side trking:
                        if exist(obj.Trkland.fx.out.trks_rh,'file') == 0
                            if strcmp(obj.projectID,'ADRC')
                                exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                    ' --seed_count=10000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                    ' --seed=' obj.Trkland.fx.in.roi_rh ' --roa=' obj.Trkland.fx.in.roa_rh_hollow ...
                                    ' --threshold_index=nqa  --fa_threshold=0.04 --fiber_count=500' ...
                                    ' --step_size=1 --turning_angle=40 --min_length=80 --max_length=250 ' ...
                                    ' --output=' obj.Trkland.fx.out.trks_rh ];
                            else
                                exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                    ' --seed_count=10000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                    ' --seed=' obj.Trkland.fx.in.roi_rh ' --roa=' obj.Trkland.fx.in.roa_rh_hollow ...
                                    ' --step_size=1 --turning_angle=40 --min_length=80 --max_length=250 ' ...
                                    ' --output=' obj.Trkland.fx.out.trks_rh ];
                            end
                            %Trying 4 times to get a trk, if not....quit
                            for dd=1:4
                                if exist(obj.Trkland.fx.out.trks_rh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            
                            wasRun=true;
                            obj.UpdateHist(obj.Trkland.fx,'trkland_fx', obj.Trkland.fx.out.trks_lh,wasRun);
                        else
                            [~, bb, cc ] = fileparts(obj.Trkland.fx.out.trks_rh);
                            fprintf(['The file ' bb cc ' exists. \n']);
                        end
                    else
                        display('QC_flag_rh found in trklnad_fx. Skipping and removing data points...')
%                         RefreshFields(obj,'fx','rh')
                    end
                else
                    display('QC_flag_bil found in trklnad_fx. Skipping and removing data points...')
%                     RefreshFields(obj,'fx','bil')
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %CLEAN UP OF THE TRACT AND EXTRACTING CENTERLINE
            for tohide=1:1
                %For left side (centerline approach):
                for tohide=1:1
                    %CLEANING UP THE STREAMLINES
                    if exist(obj.Trkland.fx.out.clean_trks_lh,'file') == 0 && exist(obj.Trkland.fx.out.trks_lh,'file') ~= 0
                        obj.Trkland.Trks.fx_raw_lh = rotrk_read(obj.Trkland.fx.out.trks_lh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'fx_lh');
                        %add Scalars
                        obj.Trkland.Trks.fx_raw_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_raw_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_raw_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_raw_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        %Trim tracts here:
                        obj.Trkland.Trks.fx_trimmed_lh = rotrk_trimmedbyTOI(obj.Trkland.Trks.fx_raw_lh, ...
                            [ {obj.Trkland.fx.in.hippo_lh} {obj.Trkland.fx.in.thalamus_lh}  ], 'fx_lh');
                        
                        %Select the HDorff centerline(first pass)
                        obj.Trkland.Trks.fx_clineinit_lh = rotrk_centerline(obj.Trkland.Trks.fx_trimmed_lh,'hausdorff');
                        %Clean up based on normality of hausdorff distance
                        obj.Trkland.Trks.fx_cleantrimmed_lh = rotrk_rm_byHDorff(obj.Trkland.Trks.fx_clineinit_lh, obj.Trkland.Trks.fx_trimmed_lh,obj.Trkland.Trks.fx_trimmed_lh);
                        %saving trimmed and trimmed_clean trks:
                        rotrk_write(obj.Trkland.Trks.fx_trimmed_lh.header,obj.Trkland.Trks.fx_trimmed_lh.sstr,obj.Trkland.fx.out.clean_trkstrimmed_lh);
                        rotrk_write(obj.Trkland.Trks.fx_cleantrimmed_lh.header,obj.Trkland.Trks.fx_cleantrimmed_lh.sstr,obj.Trkland.fx.out.clean_trks_lh);
                    end
                    
                    %HighFA:
                    if exist(obj.Trkland.fx.out.clineFA_lh_highFA,'file')==0
                        %obj.Trkland.Trks.fx_cleantrimmed_lh = rotrk_rm_bylen(obj.Trkland.Trks.fx_clineinit_lh, obj.Trkland.Trks.fx_raw_lh,obj.Trkland.Trks.fx_raw_lh);
                        %Now that the TRK is clean, lets get the high_FA and get the centerline:
                        %Pick centerline based on high_sc and FA:
                        display('executing centerline_lh for fx_highFA ... ')
                        temp_clean_trk_lh = rotrk_read(obj.Trkland.fx.out.clean_trks_lh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'fx_lh_cleantrimmed');
                        temp_clean_trk_lh = rotrk_add_sc(temp_clean_trk_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_centerline(temp_clean_trk_lh, 'high_sc','FA');
                        %Adding scalars:
                        obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        rotrk_write(obj.Trkland.Trks.fx_clinehighFA_lh.header,obj.Trkland.Trks.fx_clinehighFA_lh.sstr,obj.Trkland.fx.out.clineFA_lh_highFA );
                    end
                    
                    %HDorff:
                    if exist(obj.Trkland.fx.out.clineFA_lh_HDorff,'file')==0
                        display('executing centerline_lh for fx_hDorff ... ')
                        obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_centerline(rotrk_read(obj.Trkland.fx.out.clean_trks_lh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'fx_lh_cleantrimmed'), 'hausdorff');
                        %Add scalars:
                        obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');

                        %save trks:
                        rotrk_write(obj.Trkland.Trks.fx_clineHDorff_lh.header,obj.Trkland.Trks.fx_clineHDorff_lh.sstr,obj.Trkland.fx.out.clineFA_lh_HDorff );
                    end
                    
                    %Adding improved projectivity FA for newer centerline
                    %(and values)
                    
                    %Get data of unclean trks:
                    if  exist(obj.Trkland.fx.out.trks_lh,'file') ~= 0
                        obj.Trkland.fx.data.lh_unclean_vol = obj.Trkland.Trks.fx_raw_lh.num_uvox;
                        obj.Trkland.fx.data.lh_unclean_FA = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,4));
                        obj.Trkland.fx.data.lh_unclean_RD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,5));
                        obj.Trkland.fx.data.lh_unclean_AxD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,6));
                        obj.Trkland.fx.data.lh_unclean_MD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.lh_unclean_vol = [];
                        obj.Trkland.fx.data.lh_unclean_FA = [];
                        obj.Trkland.fx.data.lh_unclean_RD = [];
                        obj.Trkland.fx.data.lh_unclean_AxD = [];
                        obj.Trkland.fx.data.lh_unclean_MD = [];
                        
                        %Fill dependency trks to []:
                        obj.remove_trkland_fields('fx_raw_lh')
                        obj.remove_trkland_fields('fx_trimmed_lh')
                        obj.remove_trkland_fields('fx_cleantrimmed_lh')
                        obj.remove_trkland_fields('fx_clinehighFA_lh')
                        obj.remove_trkland_fields('fx_clineHDorff_lh')
                    end
                    
                    %Get data of trimmed values
                    if numel(obj.Trkland.Trks.fx_trimmed_lh.sstr) ~=0
                        obj.Trkland.fx.data.lh_trimmedclean_vol = obj.Trkland.Trks.fx_trimmed_lh.num_uvox;
                        obj.Trkland.fx.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,4));
                        obj.Trkland.fx.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,5));
                        obj.Trkland.fx.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,6));
                        obj.Trkland.fx.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.lh_trimmedclean_vol = [];
                        obj.Trkland.fx.data.lh_trimmedclean_FA = [];
                        obj.Trkland.fx.data.lh_trimmedclean_RD = [];
                        obj.Trkland.fx.data.lh_trimmedclean_AxD = [];
                        obj.Trkland.fx.data.lh_trimmedclean_MD = [];
                    end
                    
                    %Clean_fx:
                    if numel(obj.Trkland.Trks.fx_cleantrimmed_lh.sstr) ~= 0
                        obj.Trkland.fx.data.lh_clean_vol =     obj.Trkland.Trks.fx_cleantrimmed_lh.num_uvox;
                        obj.Trkland.fx.data.lh_clean_FA = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,4));
                        obj.Trkland.fx.data.lh_clean_RD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,5));
                        obj.Trkland.fx.data.lh_clean_AxD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,6));
                        obj.Trkland.fx.data.lh_clean_MD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.lh_clean_vol = [];
                        obj.Trkland.fx.data.lh_clean_FA = [];
                        obj.Trkland.fx.data.lh_clean_RD = [];
                        obj.Trkland.fx.data.lh_clean_AxD = [];
                        obj.Trkland.fx.data.lh_clean_MD = [];
                    end
                    
                    %Cline_HighFA
                    if numel(obj.Trkland.Trks.fx_clinehighFA_lh.sstr) ~= 0
                        obj.Trkland.fx.data.lh_cline_HighFA_vol = obj.Trkland.Trks.fx_clinehighFA_lh.num_uvox;
                        obj.Trkland.fx.data.lh_cline_length_highFA=obj.Trkland.Trks.fx_clinehighFA_lh.maxsstrlen;
                        obj.Trkland.fx.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,4));
                        obj.Trkland.fx.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,5));
                        obj.Trkland.fx.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,6));
                        obj.Trkland.fx.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.lh_cline_HighFA_vol = [];
                        obj.Trkland.fx.data.lh_cline_length_highFA= [];
                        obj.Trkland.fx.data.lh_cline_FA_highFA = [];
                        obj.Trkland.fx.data.lh_cline_RD_highFA = [];
                        obj.Trkland.fx.data.lh_cline_AxD_highFA = [];
                        obj.Trkland.fx.data.lh_cline_MD_highFA = [];
                    end
                    
                    %Cline_HDorff
                    if numel(obj.Trkland.Trks.fx_clineHDorff_lh.sstr) ~= 0
                        obj.Trkland.fx.data.lh_cline_HDorff_vol = obj.Trkland.Trks.fx_clineHDorff_lh.num_uvox;
                        obj.Trkland.fx.data.lh_cline_length_HDorff=obj.Trkland.Trks.fx_clineHDorff_lh.maxsstrlen;
                        obj.Trkland.fx.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,4));
                        obj.Trkland.fx.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,5));
                        obj.Trkland.fx.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,6));
                        obj.Trkland.fx.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.lh_cline_HDorff_vol = [];
                        obj.Trkland.fx.data.lh_cline_length_HDorff= [];
                        obj.Trkland.fx.data.lh_cline_FA_HDorff = [];
                        obj.Trkland.fx.data.lh_cline_RD_HDorff = [];
                        obj.Trkland.fx.data.lh_cline_AxD_HDorff = [];
                        obj.Trkland.fx.data.lh_cline_MD_HDorff = [];
                    end
                end
                
                %For right side (centerline approach):
                for tohide=1:1
                    %CLEANING UP THE STREAMLINES
                    if exist(obj.Trkland.fx.out.clean_trks_rh,'file') == 0 && exist(obj.Trkland.fx.out.trks_rh,'file') ~= 0
                        obj.Trkland.fx.out.QC = true;
                        clear obj.Trkland.Trks.fx_raw_rh obj.Trkland.Trks.fx_cleantrimmed_rh obj.Trkland.Trks.fx_clineinit_rh
                        obj.Trkland.Trks.fx_raw_rh = rotrk_read(obj.Trkland.fx.out.trks_rh,obj.sessionname ...
                            ,obj.Params.Dtifit.out.FA{end}, 'fx_rh');
                        %add Scalars
                        obj.Trkland.Trks.fx_raw_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_raw_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_raw_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_raw_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_raw_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        %Trim tracts here:
                        obj.Trkland.Trks.fx_trimmed_rh = rotrk_trimmedbyTOI(obj.Trkland.Trks.fx_raw_rh, ...
                            [ {obj.Trkland.fx.in.hippo_rh}  {obj.Trkland.fx.in.thalamus_rh}  ], 'fx_rh');
                        
                        %Select the HDorff centerline(first pass)
                        obj.Trkland.Trks.fx_clineinit_rh = rotrk_centerline(obj.Trkland.Trks.fx_trimmed_rh,'hausdorff');
                        %Clean up based on normality of hausdorff distance
                        obj.Trkland.Trks.fx_cleantrimmed_rh = rotrk_rm_byHDorff(obj.Trkland.Trks.fx_clineinit_rh, obj.Trkland.Trks.fx_trimmed_rh,obj.Trkland.Trks.fx_trimmed_rh);
                        
                        %saving trimmed and trimmed_clean trks:
                        rotrk_write(obj.Trkland.Trks.fx_trimmed_rh.header,obj.Trkland.Trks.fx_trimmed_rh.sstr,obj.Trkland.fx.out.clean_trkstrimmed_rh);
                        rotrk_write(obj.Trkland.Trks.fx_cleantrimmed_rh.header,obj.Trkland.Trks.fx_cleantrimmed_rh.sstr,obj.Trkland.fx.out.clean_trks_rh)
                    end
                    %HighFA centerline:
                    if exist(obj.Trkland.fx.out.clineFA_rh_highFA, 'file')==0
                        display('executing centerline_rh for fx_highFA ... ')
                        temp_clean_trk_rh = rotrk_read(obj.Trkland.fx.out.clean_trks_rh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'fx_rh_cleantrimmed');
                        temp_clean_trk_rh = rotrk_add_sc(temp_clean_trk_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_centerline(temp_clean_trk_rh, 'high_sc','FA');
                        obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        %save trks:
                        rotrk_write(obj.Trkland.Trks.fx_clinehighFA_rh.header,obj.Trkland.Trks.fx_clinehighFA_rh.sstr,obj.Trkland.fx.out.clineFA_rh_highFA )
                    end
                    %HDorff centerline:
                    if exist(obj.Trkland.fx.out.clineFA_rh_HDorff , 'file' ) == 0
                        display('executing centerline_rh for fx_hDorff ... ')
                        obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_centerline(rotrk_read(obj.Trkland.fx.out.clean_trks_rh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'fx_rh_cleantrimmed'), 'hausdorff');
                        %Adding scalars:
                        obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.fx_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        %save trks:
                        rotrk_write(obj.Trkland.Trks.fx_clineHDorff_rh.header,obj.Trkland.Trks.fx_clineHDorff_rh.sstr,obj.Trkland.fx.out.clineFA_rh_HDorff )
                    end
                    
                                
                    %Get data of unclean trks:
                    if  exist(obj.Trkland.fx.out.trks_rh,'file') ~= 0
                        obj.Trkland.fx.data.rh_unclean_vol = obj.Trkland.Trks.fx_raw_rh.num_uvox;
                        obj.Trkland.fx.data.rh_unclean_FA = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,4));
                        obj.Trkland.fx.data.rh_unclean_RD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,5));
                        obj.Trkland.fx.data.rh_unclean_AxD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,6));
                        obj.Trkland.fx.data.rh_unclean_MD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.rh_unclean_vol = [];
                        obj.Trkland.fx.data.rh_unclean_FA = [];
                        obj.Trkland.fx.data.rh_unclean_RD = [];
                        obj.Trkland.fx.data.rh_unclean_AxD = [];
                        obj.Trkland.fx.data.rh_unclean_MD = [];
                        
                        %Fill dependency trks to []:
                        obj.remove_trkland_fields('fx_raw_rh')
                        obj.remove_trkland_fields('fx_trimmed_rh')
                        obj.remove_trkland_fields('fx_cleantrimmed_rh')
                        obj.remove_trkland_fields('fx_clinehighFA_rh')
                        obj.remove_trkland_fields('fx_clineHDorff_rh')
                    end
                    
                    %Get data of trimmed values
                    if numel(obj.Trkland.Trks.fx_trimmed_rh.sstr) ~=0
                        obj.Trkland.fx.data.rh_trimmedclean_vol = obj.Trkland.Trks.fx_trimmed_rh.num_uvox;
                        obj.Trkland.fx.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,4));
                        obj.Trkland.fx.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,5));
                        obj.Trkland.fx.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,6));
                        obj.Trkland.fx.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.rh_trimmedclean_vol = [];
                        obj.Trkland.fx.data.rh_trimmedclean_FA = [];
                        obj.Trkland.fx.data.rh_trimmedclean_RD = [];
                        obj.Trkland.fx.data.rh_trimmedclean_AxD = [];
                        obj.Trkland.fx.data.rh_trimmedclean_MD = [];
                    end
                    
                    %Clean_fx:
                    if numel(obj.Trkland.Trks.fx_cleantrimmed_rh.sstr) ~= 0
                        obj.Trkland.fx.data.rh_clean_vol =     obj.Trkland.Trks.fx_cleantrimmed_rh.num_uvox;
                        obj.Trkland.fx.data.rh_clean_FA = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,4));
                        obj.Trkland.fx.data.rh_clean_RD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,5));
                        obj.Trkland.fx.data.rh_clean_AxD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,6));
                        obj.Trkland.fx.data.rh_clean_MD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.rh_clean_vol = [];
                        obj.Trkland.fx.data.rh_clean_FA = [];
                        obj.Trkland.fx.data.rh_clean_RD = [];
                        obj.Trkland.fx.data.rh_clean_AxD = [];
                        obj.Trkland.fx.data.rh_clean_MD = [];
                    end
                    
                    %Cline_HighFA
                    if numel(obj.Trkland.Trks.fx_clinehighFA_rh.sstr) ~= 0
                        obj.Trkland.fx.data.rh_cline_HighFA_vol = obj.Trkland.Trks.fx_clinehighFA_rh.num_uvox;
                        obj.Trkland.fx.data.rh_cline_length_highFA=obj.Trkland.Trks.fx_clinehighFA_rh.maxsstrlen;
                        obj.Trkland.fx.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,4));
                        obj.Trkland.fx.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,5));
                        obj.Trkland.fx.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,6));
                        obj.Trkland.fx.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.rh_cline_HighFA_vol = [];
                        obj.Trkland.fx.data.rh_cline_length_highFA= [];
                        obj.Trkland.fx.data.rh_cline_FA_highFA = [];
                        obj.Trkland.fx.data.rh_cline_RD_highFA = [];
                        obj.Trkland.fx.data.rh_cline_AxD_highFA = [];
                        obj.Trkland.fx.data.rh_cline_MD_highFA = [];
                    end
                    
                    %Cline_HDorff
                    if numel(obj.Trkland.Trks.fx_clineHDorff_rh.sstr) ~= 0
                        obj.Trkland.fx.data.rh_cline_HDorff_vol = obj.Trkland.Trks.fx_clineHDorff_rh.num_uvox;
                        obj.Trkland.fx.data.rh_cline_length_HDorff=obj.Trkland.Trks.fx_clineHDorff_rh.maxsstrlen;
                        obj.Trkland.fx.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,4));
                        obj.Trkland.fx.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,5));
                        obj.Trkland.fx.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,6));
                        obj.Trkland.fx.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,7));
                    else
                        obj.Trkland.fx.data.rh_cline_HDorff_vol = [];
                        obj.Trkland.fx.data.rh_cline_length_HDorff= [];
                        obj.Trkland.fx.data.rh_cline_FA_HDorff = [];
                        obj.Trkland.fx.data.rh_cline_RD_HDorff = [];
                        obj.Trkland.fx.data.rh_cline_AxD_HDorff = [];
                        obj.Trkland.fx.data.rh_cline_MD_HDorff = [];
                    end
                    
                end
            end
        end
        
        function obj = trkland_hippocing(obj)
            wasRun = false;
            fprintf('\n%s\n', 'PERFORMING TRKLAND HIPPOCAMPAL CINGULUM: TRKLAND_HIPPOCING():');
            %Create trkland directory (if doesn't exist)
            exec_cmd = [ 'mkdir -p ' obj.Trkland.root ];
            obj.RunBash(exec_cmd);
            
            %Creating root directory:
            for tohide=1:1
                outpath=obj.Trkland.root;
            end
            
            %ROIs/SEEDs PREPARATION
            for tohide=1:1
                %Dilate Hippos:
                tmp_roi_hippocing_hippo_lh = [ obj.Trkland.root 'hippocing_seed_hippoDil1_lh.nii.gz' ] ;
                tmp_roi_hippocing_hippo_rh = [ obj.Trkland.root 'hippocing_seed_hippoDil1_rh.nii.gz' ] ;
                obj.Trkland.hippocing.in.seed_hippo_lh = [ obj.Trkland.root 'hippocing_seed_hippoDil2_lh.nii.gz' ] ;
                obj.Trkland.hippocing.in.seed_hippo_rh = [ obj.Trkland.root 'hippocing_seed_hippoDil2_rh.nii.gz' ] ;
                %left:
                if exist(tmp_roi_hippocing_hippo_lh , 'file') == 0
                    fprintf('\n Dilating Hippocampus_lh...')
                    exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.hippo_lh ...
                        ' -dilM ' tmp_roi_hippocing_hippo_lh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                if exist(obj.Trkland.hippocing.in.seed_hippo_lh , 'file') == 0
                    fprintf('\n 2nd dilation Hippocampus_lh...')
                    exec_cmd = ['fslmaths ' tmp_roi_hippocing_hippo_lh ...
                        ' -dilM ' obj.Trkland.hippocing.in.seed_hippo_lh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %right:
                if exist( tmp_roi_hippocing_hippo_rh , 'file') == 0
                    fprintf('\n Dilating Hippocampus_rh...')
                    exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.hippo_rh ...
                        ' -dilM '  tmp_roi_hippocing_hippo_rh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                if exist(obj.Trkland.hippocing.in.seed_hippo_rh , 'file') == 0
                    fprintf('\n 2nd dilation Hippocampus_rh...')
                    exec_cmd = ['fslmaths ' tmp_roi_hippocing_hippo_rh  ...
                        ' -dilM ' obj.Trkland.hippocing.in.seed_hippo_rh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                
                %Dilate Postcingulates:
                if strcmp(obj.projectID,'HAB')
                    tmp_roi_hippocing_postcing_lh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil1_lh.nii.gz' ] ;
                    tmp_roi_hippocing_postcing_rh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil1_rh.nii.gz' ] ;
                    obj.Trkland.hippocing.in.roi_postcing_lh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil2_lh.nii.gz' ] ;
                    obj.Trkland.hippocing.in.roi_postcing_rh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil2_rh.nii.gz' ] ;
                    
                    %left:
                    if exist(tmp_roi_hippocing_postcing_lh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_lh...')
                        exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.postcing_lh ...
                            ' -dilM ' tmp_roi_hippocing_postcing_lh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                    if exist(obj.Trkland.hippocing.in.roi_postcing_lh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_lh...')
                        exec_cmd = ['fslmaths '  tmp_roi_hippocing_postcing_lh ...
                            ' -dilM ' obj.Trkland.hippocing.in.roi_postcing_lh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                    
                    %right:
                    if exist(tmp_roi_hippocing_postcing_rh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_lh...')
                        exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.postcing_lh ...
                            ' -dilM ' tmp_roi_hippocing_postcing_rh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                    if exist(obj.Trkland.hippocing.in.roi_postcing_rh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_rh...')
                        exec_cmd = ['fslmaths ' tmp_roi_hippocing_postcing_rh ...
                            ' -dilM ' obj.Trkland.hippocing.in.roi_postcing_rh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                else
                    
                    obj.Trkland.hippocing.in.roi_postcing_lh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil1_lh.nii.gz' ] ;
                    obj.Trkland.hippocing.in.roi_postcing_rh = [ obj.Trkland.root 'hippocing_roi_postcingulateDil1_rh.nii.gz' ] ;
                    
                    %left:
                    if exist(obj.Trkland.hippocing.in.roi_postcing_lh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_lh...')
                        exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.postcing_lh ...
                            ' -dilM ' obj.Trkland.hippocing.in.roi_postcing_lh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                    
                    %right:
                    if exist(obj.Trkland.hippocing.in.roi_postcing_rh  , 'file') == 0
                        fprintf('\n Dilating Posterior cingulate_rh...')
                        exec_cmd = ['fslmaths ' obj.Trkland.hippocing.in.postcing_rh ...
                            ' -dilM ' obj.Trkland.hippocing.in.roi_postcing_rh   ];
                        obj.RunBash(exec_cmd);
                        fprintf('...done \n');
                    end
                end
            end
            %INIT CLEANUP OF THE TRACTS
            for tohide=1:1
                obj.Trkland.hippocing.out.clean_trkstrimmed_lh = [ obj.Trkland.root  'trkk_hippocing_trimmed_lh.trk.gz'];
                obj.Trkland.hippocing.out.clean_trkstrimmed_rh = [ obj.Trkland.root  'trkk_hippocing_trimmed_rh.trk.gz'];
                
                obj.Trkland.hippocing.out.clean_trks_lh = [ obj.Trkland.root  'trkk_hippocing_trimmedclean_lh.trk.gz'];
                obj.Trkland.hippocing.out.clean_trks_rh = [ obj.Trkland.root  'trkk_hippocing_trimmedclean_rh.trk.gz'];
                
                obj.Trkland.hippocing.out.clineFA_lh_highFA = [ obj.Trkland.root  'cline_hippocing_highFA_lh.trk.gz'];
                obj.Trkland.hippocing.out.clineFA_lh_HDorff = [ obj.Trkland.root  'cline_hippocing_HDorff_lh.trk.gz'];
                obj.Trkland.hippocing.out.clineFA_rh_highFA = [ obj.Trkland.root  'cline_hippocing_highFA_rh.trk.gz'];
                obj.Trkland.hippocing.out.clineFA_rh_HDorff = [ obj.Trkland.root  'cline_hippocing_HDorff_rh.trk.gz'];
                
                obj.Trkland.hippocing.QCfile_lh = [outpath 'QC_hcing_lh.flag'] ;
                obj.Trkland.hippocing.QCfile_rh = [outpath 'QC_hcing_rh.flag'] ;
                obj.Trkland.hippocing.QCfile_bil = [outpath 'QC_hcing_bil.flag'] ;
            end
            
            %TRACKING STARS HERE:
            for tohide=1:1
                obj.Trkland.hippocing.out.trk_lh = [ obj.Trkland.root  'trkk_hippocing_lh.trk.gz'];
                obj.Trkland.hippocing.out.trk_rh = [ obj.Trkland.root  'trkk_hippocing_rh.trk.gz'];
                 if exist(obj.Trkland.hippocing.QCfile_bil, 'file') == 0
                     %Left side trking:
                    if exist(obj.Trkland.hippocing.QCfile_lh, 'file') == 0
                        if exist(obj.Trkland.hippocing.out.trk_lh,'file') == 0
                            exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                ' --seed=' obj.Trkland.hippocing.in.seed_hippo_lh ' --roi=' obj.Trkland.hippocing.in.roi_postcing_lh ...
                                ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
                                ' --output=' obj.Trkland.hippocing.out.trk_lh ];
                            for dd=1:4 %trying 4 times to get a trk. If not, quit!
                                if exist(obj.Trkland.hippocing.out.trk_lh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            wasRun=true;
                            obj.UpdateHist(obj.Trkland.hippocing,'trkland_hippocing', obj.Trkland.hippocing.out.trk_lh,wasRun);
                        end
                    else
                        display('QC_flag_lh found in trkland_hippocing. Skipping tracking...')
                    end
                    
                    %Right side trking:
                    if exist(obj.Trkland.hippocing.QCfile_rh, 'file') == 0
                        if exist(obj.Trkland.hippocing.out.trk_rh,'file') == 0
                            exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                ' --seed=' obj.Trkland.hippocing.in.seed_hippo_rh ' --roi=' obj.Trkland.hippocing.in.roi_postcing_rh ...
                                ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
                                ' --output=' obj.Trkland.hippocing.out.trk_rh ];
                            for dd=1:4 %trying 4 times to get a trk. If not, quit!
                                if exist(obj.Trkland.hippocing.out.trk_rh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            wasRun=true;
                            obj.UpdateHist(obj.Trkland.hippocing,'trkland_hippocing', obj.Trkland.hippocing.out.trk_rh,wasRun);
                        end
                    else
                        display('QC_flag_rh found in trkland_hippocing. Skipping tracking...')
                    end
                else
                    display('QC_flag_bil found in trkland_hippocing. Skiping tracking...')
                end
            end
            
            %LEFT SIDE:
            for tohide=1:1
                if exist(obj.Trkland.hippocing.out.clean_trks_lh ,'file') == 0 && exist(obj.Trkland.hippocing.out.trk_lh,'file') ~= 0
                    obj.Trkland.Trks.raw_hippocing_lh = rotrk_read(obj.Trkland.hippocing.out.trk_lh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'hippocing_lh');
                    %add Scalars:
                    obj.Trkland.Trks.raw_hippocing_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.raw_hippocing_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.raw_hippocing_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.raw_hippocing_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %Trim tracts here:
                    obj.Trkland.Trks.hippocing_trimmed_lh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_hippocing_lh, ...
                        [  {obj.Trkland.hippocing.in.hippo_lh}  {obj.Trkland.hippocing.in.roi_postcing_lh}  ], 'postcing_lh');
                    %Select the HDorff centerline(first pass):
                    obj.Trkland.Trks.hippocing_clineinit_lh= rotrk_centerline(obj.Trkland.Trks.hippocing_trimmed_lh,'hausdorff');
                    %Clean up based on normality of hausdorff distance:
                    obj.Trkland.Trks.hippocing_cleantrimmed_lh = rotrk_rm_byHDorff(obj.Trkland.Trks.hippocing_clineinit_lh, obj.Trkland.Trks.hippocing_trimmed_lh,obj.Trkland.Trks.hippocing_trimmed_lh);
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_trimmed_lh.header,obj.Trkland.Trks.hippocing_trimmed_lh.sstr,obj.Trkland.hippocing.out.clean_trkstrimmed_lh);
                    rotrk_write(obj.Trkland.Trks.hippocing_cleantrimmed_lh.header,obj.Trkland.Trks.hippocing_cleantrimmed_lh.sstr,obj.Trkland.hippocing.out.clean_trks_lh );
                end
                
                %HighFA_centerline: 
                if exist(obj.Trkland.hippocing.out.clineFA_lh_highFA,'file')==0
                    display('Executing centerline_lh for hippocing_highFA...')
                    temp_clean_trk_lh=rotrk_read(obj.Trkland.hippocing.out.clean_trks_lh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'hippocing_lh_cleantrimmed');
                    temp_clean_trk_lh=rotrk_add_sc(temp_clean_trk_lh,obj.Params.Dtifit.out.FA{end},'FA');
                    %Adding scalars:
                    obj.Trkland.Trks.hippocing_clinehighFA_lh = rotrk_centerline(temp_clean_trk_lh, 'high_sc','FA');
                    obj.Trkland.Trks.hippocing_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.hippocing_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.hippocing_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.hippocing_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_clinehighFA_lh.header,obj.Trkland.Trks.hippocing_clinehighFA_lh.sstr,obj.Trkland.hippocing.out.clineFA_lh_highFA );
                end
                %HDorff_centerline:
                if exist(obj.Trkland.hippocing.out.clineFA_lh_HDorff,'file')==0
                    obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_centerline(rotrk_read(obj.Trkland.hippocing.out.clean_trks_lh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'hippocing_lh_cleantrimmed'), 'hausdorff');
                    %Adding scalars:
                    obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_clineHDorff_lh.header,obj.Trkland.Trks.hippocing_clineHDorff_lh.sstr,obj.Trkland.hippocing.out.clineFA_lh_HDorff);
                end
                
                %GETTING DATA NOW:
                %Unclean tracts:
                if exist(obj.Trkland.hippocing.out.trk_lh,'file') ~= 0
                    obj.Trkland.hippocing.data.lh_unclean_vol = obj.Trkland.Trks.raw_hippocing_lh.num_uvox;
                    obj.Trkland.hippocing.data.lh_unclean_FA = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.lh_unclean_RD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.lh_unclean_AxD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.lh_unclean_MD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.lh_unclean_vol =  []; 
                    obj.Trkland.hippocing.data.lh_unclean_FA = [];
                    obj.Trkland.hippocing.data.lh_unclean_RD = [];
                    obj.Trkland.hippocing.data.lh_unclean_AxD = [];
                    obj.Trkland.hippocing.data.lh_unclean_MD = [];
                    
                    %Fill dependenc trks to[]:
                    obj.remove_trkland_fields('raw_hippocing_lh')
                    obj.remove_trkland_fields('hippocing_trimmed_lh')
                    obj.remove_trkland_fields('hippocing_cleantrimmed_lh')
                    obj.remove_trkland_fields('hippocing_clinehighFA_lh')
                    obj.remove_trkland_fields('hippocing_clineHDorff_lh')
                 end
                
                
                %Trimmed_clean_hippocing:
                if numel(obj.Trkland.Trks.hippocing_trimmed_lh.sstr) ~= 0
                    obj.Trkland.hippocing.data.lh_trimmedclean_vol = obj.Trkland.Trks.hippocing_trimmed_lh.num_uvox;
                    obj.Trkland.hippocing.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.lh_trimmedclean_vol = [];
                    obj.Trkland.hippocing.data.lh_trimmedclean_FA = [];
                    obj.Trkland.hippocing.data.lh_trimmedclean_RD = [];
                    obj.Trkland.hippocing.data.lh_trimmedclean_AxD = [];
                    obj.Trkland.hippocing.data.lh_trimmedclean_MD = [];
                end
                
                %Clean_hippocing:
                if numel(obj.Trkland.Trks.hippocing_cleantrimmed_lh.sstr) ~= 0
                    obj.Trkland.hippocing.data.lh_clean_vol = obj.Trkland.Trks.hippocing_cleantrimmed_lh.num_uvox;
                    obj.Trkland.hippocing.data.lh_clean_FA = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.lh_clean_RD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.lh_clean_AxD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.lh_clean_MD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,7));
                else
                    
                    obj.Trkland.hippocing.data.lh_clean_vol = [];
                    obj.Trkland.hippocing.data.lh_clean_FA = [];
                    obj.Trkland.hippocing.data.lh_clean_RD = [];
                    obj.Trkland.hippocing.data.lh_clean_AxD = [];
                    obj.Trkland.hippocing.data.lh_clean_MD = [];
                end
                
                %Cline_HighFA
                if numel(obj.Trkland.Trks.hippocing_clinehighFA_lh.sstr) ~= 0
                    obj.Trkland.hippocing.data.lh_cline_HighFA_vol = obj.Trkland.Trks.hippocing_clinehighFA_lh.num_uvox;
                    obj.Trkland.hippocing.data.lh_cline_length_highFA=obj.Trkland.Trks.hippocing_clinehighFA_lh.maxsstrlen;
                    obj.Trkland.hippocing.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.lh_cline_HighFA_vol = [];
                    obj.Trkland.hippocing.data.lh_cline_length_highFA= [];
                    obj.Trkland.hippocing.data.lh_cline_FA_highFA = [];
                    obj.Trkland.hippocing.data.lh_cline_RD_highFA = [];
                    obj.Trkland.hippocing.data.lh_cline_AxD_highFA = [];
                    obj.Trkland.hippocing.data.lh_cline_MD_highFA = [];
                end
                
                %Cline_HDorff
                if numel(obj.Trkland.Trks.hippocing_clineHDorff_lh.sstr) ~= 0
                    obj.Trkland.hippocing.data.lh_cline_HDorff_vol = obj.Trkland.Trks.hippocing_clineHDorff_lh.num_uvox;
                    obj.Trkland.hippocing.data.lh_cline_length_HDorff=obj.Trkland.Trks.hippocing_clineHDorff_lh.maxsstrlen;
                    obj.Trkland.hippocing.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.lh_cline_HDorff_vol = [];
                    obj.Trkland.hippocing.data.lh_cline_length_HDorff = [];
                    obj.Trkland.hippocing.data.lh_cline_FA_HDorff = [];
                    obj.Trkland.hippocing.data.lh_cline_RD_HDorff = [];
                    obj.Trkland.hippocing.data.lh_cline_AxD_HDorff = [];
                    obj.Trkland.hippocing.data.lh_cline_MD_HDorff = [];
                    
                end
            end
            
            %RIGHT SIDE:
            for tohide=1:1
                if exist(obj.Trkland.hippocing.out.clean_trks_rh ,'file') == 0 && exist(obj.Trkland.hippocing.out.trk_rh,'file') ~= 0
                    
                    clear obj.Trkland.Trks.raw_hippocing_rh obj.Trkland.Trks.hippocing_cleantrimmed_rh obj.Trkland.Trks.hippocing_clineinit_rh
                    obj.Trkland.Trks.raw_hippocing_rh = rotrk_read(obj.Trkland.hippocing.out.trk_rh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'hippocing_rh');
                    %add Scalars
                    obj.Trkland.Trks.raw_hippocing_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.raw_hippocing_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.raw_hippocing_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.raw_hippocing_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_hippocing_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %Trim tracts here:
                    if exist(obj.Trkland.hippocing.out.clean_trkstrimmed_rh,'file') == 0
                        obj.Trkland.Trks.hippocing_trimmed_rh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_hippocing_rh, ...
                            [ {obj.Trkland.hippocing.in.hippo_rh}   {obj.Trkland.hippocing.in.roi_postcing_rh}   ], 'postcing_rh');
                    else
                        obj.Trkland.Trks.hippocing_trimmed_rh =    rotrk_read(obj.Trkland.hippocing.out.clean_trkstrimmed_rh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'hippocing_lh' );
                        obj.Trkland.Trks.hippocing_trimmed_rh  = rotrk_add_sc(obj.Trkland.Trks.hippocing_trimmed_rh  ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.hippocing_trimmed_rh  = rotrk_add_sc(obj.Trkland.Trks.hippocing_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.hippocing_trimmed_rh = rotrk_add_sc(obj.Trkland.Trks.hippocing_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.hippocing_trimmed_rh  = rotrk_add_sc( obj.Trkland.Trks.hippocing_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        
                    end
                    %
                    %Select the HDorff centerline(first pass)
                    obj.Trkland.Trks.hippocing_clineinit_rh= rotrk_centerline(obj.Trkland.Trks.hippocing_trimmed_rh,'hausdorff');
                    %Clean up based on normality of hausdorff distance
                    obj.Trkland.Trks.hippocing_cleantrimmed_rh = rotrk_rm_byHDorff(obj.Trkland.Trks.hippocing_clineinit_rh, obj.Trkland.Trks.hippocing_trimmed_rh,obj.Trkland.Trks.hippocing_trimmed_rh);
                    %%%%obj.Trkland.Trks.hippocing_cleantrimmed_rh = rotrk_rm_bylen(temp_fx_rh_cline, temp_fx_rh,temp_fx_rh);
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_cleantrimmed_rh.header,obj.Trkland.Trks.hippocing_cleantrimmed_rh.sstr,obj.Trkland.hippocing.out.clean_trks_rh )
                    rotrk_write(obj.Trkland.Trks.hippocing_trimmed_rh.header,obj.Trkland.Trks.hippocing_trimmed_rh.sstr,obj.Trkland.hippocing.out.clean_trkstrimmed_rh);
                    
                    wasRun=true;
                    obj.UpdateHist(obj.Trkland.hippocing,'trkland_hippocing', obj.Trkland.hippocing.out.clineFA_rh_highFA,wasRun);
                    
                end
                 
                %Now that the TRK is clean, lets get the high_FA and get the centerline:
                %Pick centerline based on high_sc and FA:
                %HighFA
                if exist(obj.Trkland.hippocing.out.clineFA_rh_highFA ,'file') == 0
                    display('Executing centerline_rh for hippocing_highFA...');
                    temp_clean_trk_rh=rotrk_read(obj.Trkland.hippocing.out.clean_trks_rh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'hippocing_rh_cleantrimmed');
                    temp_clean_trk_rh=rotrk_add_sc(temp_clean_trk_rh,obj.Params.Dtifit.out.FA{end},'FA');
                    obj.Trkland.Trks.hippocing_clinehighFA_rh= rotrk_centerline(temp_clean_trk_rh, 'high_sc','FA');
                    %Adding scalars:
                    obj.Trkland.Trks.hippocing_clinehighFA_rh = rotrk_centerline(temp_clean_trk_rh, 'high_sc','FA');
                    obj.Trkland.Trks.hippocing_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.hippocing_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.hippocing_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.hippocing_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_clinehighFA_rh.header,obj.Trkland.Trks.hippocing_clinehighFA_rh.sstr,obj.Trkland.hippocing.out.clineFA_rh_highFA )
                end
                %Hdorff
                if exist(obj.Trkland.hippocing.out.clineFA_rh_HDorff, 'file') == 0
                    obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_centerline(rotrk_read(obj.Trkland.hippocing.out.clean_trks_rh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'hippocing_rh_cleantrimmed'), 'hausdorff');
                    %Adding scalars:
                    obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.hippocing_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.hippocing_clineHDorff_rh.header,obj.Trkland.Trks.hippocing_clineHDorff_rh.sstr,obj.Trkland.hippocing.out.clineFA_rh_HDorff)
                end
                
                %Get info about raw data:
                if exist(obj.Trkland.hippocing.out.trk_rh,'file') == 0 %Either failed
                    obj.Trkland.Trks.raw_hippocing_rh.sstr = [];
                    obj.Trkland.Trks.raw_hippocing_rh.unique_voxels = [];
                    
                end
                %Get volume data of unclean/cleaned tracts:
                %Unclean tracts:
                if exist(obj.Trkland.hippocing.out.trk_rh,'file') ~= 0
                    obj.Trkland.hippocing.data.rh_unclean_vol = obj.Trkland.Trks.raw_hippocing_rh.num_uvox;
                    obj.Trkland.hippocing.data.rh_unclean_FA = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.rh_unclean_RD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.rh_unclean_AxD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.rh_unclean_MD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.rh_unclean_vol = [];
                    obj.Trkland.hippocing.data.rh_unclean_FA = [];
                    obj.Trkland.hippocing.data.rh_unclean_RD = [];
                    obj.Trkland.hippocing.data.rh_unclean_AxD = [];
                    obj.Trkland.hippocing.data.rh_unclean_MD = [];
                    
                    %Fill dependenc trks to[]:
                    obj.remove_trkland_fields('raw_hippocing_rh')
                    obj.remove_trkland_fields('hippocing_trimmed_rh')
                    obj.remove_trkland_fields('hippocing_cleantrimmed_rh')
                    obj.remove_trkland_fields('hippocing_clinehighFA_rh')
                    obj.remove_trkland_fields('hippocing_clineHDorff_rh')
                end
                
                %Trimmed_clean_hippocing:
                if numel(obj.Trkland.Trks.hippocing_trimmed_rh.sstr) ~= 0
                    obj.Trkland.hippocing.data.rh_trimmedclean_vol = obj.Trkland.Trks.hippocing_trimmed_rh.num_uvox;
                    obj.Trkland.hippocing.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.rh_trimmedclean_vol = [];
                    obj.Trkland.hippocing.data.rh_trimmedclean_FA = [];
                    obj.Trkland.hippocing.data.rh_trimmedclean_RD = [];
                    obj.Trkland.hippocing.data.rh_trimmedclean_AxD = [];
                    obj.Trkland.hippocing.data.rh_trimmedclean_MD = [];
                end
                
                %Clean_hippocing:
                if numel(obj.Trkland.Trks.hippocing_cleantrimmed_rh.sstr) ~= 0
                    obj.Trkland.hippocing.data.rh_clean_vol = obj.Trkland.Trks.hippocing_cleantrimmed_rh.num_uvox;
                    obj.Trkland.hippocing.data.rh_clean_FA = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.rh_clean_RD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.rh_clean_AxD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.rh_clean_MD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,7));
                else
                    
                    obj.Trkland.hippocing.data.rh_clean_vol = [];
                    obj.Trkland.hippocing.data.rh_clean_FA = [];
                    obj.Trkland.hippocing.data.rh_clean_RD = [];
                    obj.Trkland.hippocing.data.rh_clean_AxD = [];
                    obj.Trkland.hippocing.data.rh_clean_MD = [];
                end
                
                %Cline_HighFA
                if numel(obj.Trkland.Trks.hippocing_clinehighFA_rh.sstr) ~= 0
                    obj.Trkland.hippocing.data.rh_cline_HighFA_vol = obj.Trkland.Trks.hippocing_clinehighFA_rh.num_uvox;
                    obj.Trkland.hippocing.data.rh_cline_length_highFA=obj.Trkland.Trks.hippocing_clinehighFA_rh.maxsstrlen;
                    obj.Trkland.hippocing.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.rh_cline_HighFA_vol = [];
                    obj.Trkland.hippocing.data.rh_cline_length_highFA= [];
                    obj.Trkland.hippocing.data.rh_cline_FA_highFA = [];
                    obj.Trkland.hippocing.data.rh_cline_RD_highFA = [];
                    obj.Trkland.hippocing.data.rh_cline_AxD_highFA = [];
                    obj.Trkland.hippocing.data.rh_cline_MD_highFA = [];
                end
                
                %Cline_HDorff
                if numel(obj.Trkland.Trks.hippocing_clineHDorff_rh.sstr) ~= 0
                    obj.Trkland.hippocing.data.rh_cline_HDorff_vol = obj.Trkland.Trks.hippocing_clineHDorff_rh.num_uvox;
                    obj.Trkland.hippocing.data.rh_cline_length_HDorff=obj.Trkland.Trks.hippocing_clineHDorff_rh.maxsstrlen;
                    obj.Trkland.hippocing.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,4));
                    obj.Trkland.hippocing.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,5));
                    obj.Trkland.hippocing.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,6));
                    obj.Trkland.hippocing.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,7));
                else
                    obj.Trkland.hippocing.data.rh_cline_HDorff_vol = [];
                    obj.Trkland.hippocing.data.rh_cline_length_HDorff = [];
                    obj.Trkland.hippocing.data.rh_cline_FA_HDorff = [];
                    obj.Trkland.hippocing.data.rh_cline_RD_HDorff = [];
                    obj.Trkland.hippocing.data.rh_cline_AxD_HDorff = [];
                    obj.Trkland.hippocing.data.rh_cline_MD_HDorff = [];
                end
            end
        end
        
        function obj = trkland_cingulum(obj)
            wasRun = false;
            fprintf('\n%s\n', 'PERFORMING TRKLAND CINGULUM: TRKLAND_CINGULUM():');
            
            %Create trkland directory (if doesn't exist)
            exec_cmd = [ 'mkdir -p ' obj.Trkland.root ];
            obj.RunBash(exec_cmd);
            outpath=obj.Trkland.root;
            
            %ROIs/SEEDs PREPARATION
            tmp_roi_antroscing_lh = [ obj.Trkland.root 'cingulum_roi_antrostralcingulateDil1_lh.nii.gz' ] ;
            obj.Trkland.cingulum.in.roi_antroscing_lh = [ obj.Trkland.root 'cingulum_roi_antrostralcingulateDil2_lh.nii.gz' ] ;
            
            tmp_roi_antroscing_rh = [ obj.Trkland.root 'cingulum_roi_antrostralcingulateDil1_rh.nii.gz' ] ;
            obj.Trkland.cingulum.in.roi_antroscing_rh = [ obj.Trkland.root 'cingulum_roi_antrostralcingulateDil2_rh.nii.gz' ] ;
            
            obj.Trkland.cingulum.in.seed_postcing_lh = [ obj.Trkland.root 'cingulum_seed_postcingulate_lh.nii.gz' ] ;
            obj.Trkland.cingulum.in.seed_postcing_rh = [ obj.Trkland.root 'cingulum_seed_postcingulate_rh.nii.gz' ] ;
            for tohide=1:1
                %left:
                if exist(obj.Trkland.cingulum.in.roi_antroscing_lh , 'file') == 0
                    fprintf('\n Working on anterior cingulate_lh...')
                    exec_cmd = ['fslmaths '  obj.Trkland.cingulum.in.rostantcing_lh ...
                        ' -dilM ' tmp_roi_antroscing_lh  ];
                    obj.RunBash(exec_cmd);
                    exec_cmd = ['fslmaths '  tmp_roi_antroscing_lh ...
                        ' -dilM ' obj.Trkland.cingulum.in.roi_antroscing_lh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                if exist(obj.Trkland.cingulum.in.seed_postcing_lh , 'file') == 0
                    fprintf('\n Working on posterior cingulate_lh...')
                    exec_cmd = ['cp ' obj.Trkland.cingulum.in.postcing_lh ...
                        ' ' obj.Trkland.cingulum.in.seed_postcing_lh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                %right:
                if exist(obj.Trkland.cingulum.in.roi_antroscing_rh , 'file') == 0
                    fprintf('\n Working on anterior cingulate_rh...')
                    exec_cmd = ['fslmaths '  obj.Trkland.cingulum.in.rostantcing_rh ...
                        ' -dilM ' tmp_roi_antroscing_rh  ];
                    obj.RunBash(exec_cmd);
                    exec_cmd = ['fslmaths '  tmp_roi_antroscing_rh ...
                        ' -dilM ' obj.Trkland.cingulum.in.roi_antroscing_rh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
                if exist(obj.Trkland.cingulum.in.seed_postcing_rh , 'file') == 0
                    fprintf('\n Working on posterior cingulate_rh...')
                    exec_cmd = ['cp ' obj.Trkland.cingulum.in.postcing_rh ...
                        ' ' obj.Trkland.cingulum.in.seed_postcing_rh  ];
                    obj.RunBash(exec_cmd);
                    fprintf('...done \n');
                end
            end
            
            %INIT CLEANUP VARIABLES
            for tohide=1:1
                obj.Trkland.cingulum.out.clean_trkstrimmed_lh = [ obj.Trkland.root  'trkk_cingulum_trimmed_lh.trk.gz'];
                obj.Trkland.cingulum.out.clean_trkstrimmed_rh = [ obj.Trkland.root  'trkk_cingulum_trimmed_rh.trk.gz'];
                
                obj.Trkland.cingulum.out.clean_trks_lh = [ obj.Trkland.root  'trkk_cingulum_trimmedclean_lh.trk.gz'];
                obj.Trkland.cingulum.out.clean_trks_rh = [ obj.Trkland.root  'trkk_cingulum_trimmedclean_rh.trk.gz'];
                
                obj.Trkland.cingulum.out.clineFA_lh_highFA = [ obj.Trkland.root  'cline_cingulum_highFAlh.trk.gz'];
                obj.Trkland.cingulum.out.clineFA_lh_HDorff = [ obj.Trkland.root  'cline_cingulum_HDorfflh.trk.gz'];
                obj.Trkland.cingulum.out.clineFA_rh_highFA = [ obj.Trkland.root  'cline_cingulum_highFArh.trk.gz'];
                obj.Trkland.cingulum.out.clineFA_rh_HDorff = [ obj.Trkland.root  'cline_cingulum_HDorffrh.trk.gz'];
                
                obj.Trkland.cingulum.QCfile_lh = [outpath 'QC_cingulum_lh.flag'] ;
                obj.Trkland.cingulum.QCfile_rh = [outpath 'QC_cingulum_rh.flag'] ;
                obj.Trkland.cingulum.QCfile_bil = [outpath 'QC_cingulum_bil.flag'] ;
                
            end
            
            %TRACKING STARS HERE:
            for tohide=1:1
                obj.Trkland.cingulum.out.trk_lh = [ obj.Trkland.root  'trkk_cingulum_lh.trk.gz'];
                obj.Trkland.cingulum.out.trk_rh = [ obj.Trkland.root  'trkk_cingulum_rh.trk.gz'];
                if exist(obj.Trkland.cingulum.QCfile_bil, 'file') == 0
                    %Left side trking:
                    if exist(obj.Trkland.cingulum.QCfile_lh, 'file') == 0
                        if exist(obj.Trkland.cingulum.out.trk_lh,'file') == 0
                            exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                ' --seed=' obj.Trkland.cingulum.in.seed_postcing_lh ' --roi=' obj.Trkland.cingulum.in.roi_antroscing_lh ...
                                ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
                                ' --output=' obj.Trkland.cingulum.out.trk_lh ];
                            for dd=1:4 %trying 4 times to get a trk. If not, quit!
                                if exist(obj.Trkland.cingulum.out.trk_lh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            wasRun=true;
                            obj.UpdateHist(obj.Trkland.cingulum,'trkland_cingulum', obj.Trkland.cingulum.out.trk_lh,wasRun);
                        end
                    else
                        display('QC_flag_lh found in trkland_cingulum. Skipping and removing data points...')
                   end
                    
                    %Right side trking:
                    if exist(obj.Trkland.cingulum.QCfile_rh, 'file') == 0
                        
                        if exist(obj.Trkland.cingulum.out.trk_rh,'file') == 0
                            exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
                                ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
                                ' --seed=' obj.Trkland.cingulum.in.seed_postcing_rh ' --roi=' obj.Trkland.cingulum.in.roi_antroscing_rh ...
                                ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
                                ' --output=' obj.Trkland.cingulum.out.trk_rh ];
                            
                            for dd=1:4 %trying 4 times to get a trk. If not, quit!
                                if exist(obj.Trkland.cingulum.out.trk_rh,'file') == 0
                                    obj.RunBash(exec_cmd,144);
                                end
                            end
                            wasRun=true;
                             obj.UpdateHist(obj.Trkland.cingulum,'trkland_cingulum', obj.Trkland.cingulum.out.trk_rh,wasRun);
                        end
                    else
                        display('QC_flag_rh found in trkland_cingulum. Skipping and removing data points...')
                   end
                else
                    display('QC_flag_bil found in trkland_cingulum. Skipping and removing data points...')
               end
            end
            
            %LEFT SIDE:
            for tohide=1:1
                if exist(obj.Trkland.cingulum.out.clean_trks_lh ,'file') == 0 && exist(obj.Trkland.cingulum.out.trk_lh,'file') ~= 0
                        obj.Trkland.Trks.raw_cingulum_lh = rotrk_read(obj.Trkland.cingulum.out.trk_lh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'cingulum_lh_cleantrimmed');
                        %add Scalars:
                        obj.Trkland.Trks.raw_cingulum_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.raw_cingulum_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.raw_cingulum_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.raw_cingulum_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        
                        %Trim tracts here (checking if trimming occurs, mainly for manual edits):
                        if exist(obj.Trkland.cingulum.out.clean_trkstrimmed_lh,'file') == 0
                            obj.Trkland.Trks.cingulum_trimmed_lh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_cingulum_lh, ...
                                [ {obj.Trkland.cingulum.in.seed_postcing_lh}  {obj.Trkland.cingulum.in.rostantcing_lh}  ], 'cingulum_lh');
                        else
                            obj.Trkland.Trks.cingulum_trimmed_lh = rotrk_read(obj.Trkland.cingulum.out.clean_trkstrimmed_lh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'cingulum_lh_cleantrimmed' );
                            obj.Trkland.Trks.cingulum_trimmed_lh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_lh  ,obj.Params.Dtifit.out.FA{end} , 'FA');
                            obj.Trkland.Trks.cingulum_trimmed_lh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                            obj.Trkland.Trks.cingulum_trimmed_lh = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_lh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                            obj.Trkland.Trks.cingulum_trimmed_lh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_lh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                        end
                        
                        
                        %Select the HDorff centerline(first pass):
                        obj.Trkland.Trks.cingulum_clineinit_lh= rotrk_centerline(obj.Trkland.Trks.cingulum_trimmed_lh,'hausdorff');
                        %Clean up based on normality of hausdorff distance:
                        obj.Trkland.Trks.cingulum_cleantrimmed_lh = rotrk_rm_byHDorff(obj.Trkland.Trks.cingulum_clineinit_lh, obj.Trkland.Trks.cingulum_trimmed_lh,obj.Trkland.Trks.cingulum_trimmed_lh);
                        %save trks:
                        rotrk_write(obj.Trkland.Trks.cingulum_trimmed_lh.header,obj.Trkland.Trks.cingulum_trimmed_lh.sstr,obj.Trkland.cingulum.out.clean_trkstrimmed_lh);
                        rotrk_write(obj.Trkland.Trks.cingulum_cleantrimmed_lh.header,obj.Trkland.Trks.cingulum_cleantrimmed_lh.sstr,obj.Trkland.cingulum.out.clean_trks_lh );
                        
                        wasRun=true;
                        obj.UpdateHist(obj.Trkland.cingulum,'trkland_cingulum', obj.Trkland.cingulum.out.clineFA_lh_highFA,wasRun);
                end
                %Now that the TRK is clean, lets get the high_FA and get the centerline:
                %Pick centerline based on high_sc and FA:
                
                %centerline_HighFA:
                if exist(obj.Trkland.cingulum.out.clineFA_lh_highFA,'file') == 0
                    display('Executing centerline_lh for cingulum_highFA...');
                    temp_clean_trk_lh=rotrk_read(obj.Trkland.cingulum.out.clean_trks_lh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'cingulum_lh_cleantrimmed');
                    temp_clean_trk_lh=rotrk_add_sc(temp_clean_trk_lh,obj.Params.Dtifit.out.FA{end},'FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_lh= rotrk_centerline(temp_clean_trk_lh, 'high_sc','FA');
                    %Adding scalars:
                    obj.Trkland.Trks.cingulum_clinehighFA_lh = rotrk_centerline(temp_clean_trk_lh, 'high_sc','FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.cingulum_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.cingulum_clinehighFA_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.cingulum_clinehighFA_lh.header,obj.Trkland.Trks.cingulum_clinehighFA_lh.sstr,obj.Trkland.cingulum.out.clineFA_lh_highFA );
                end
                %centerline_HDorff:
                if exist(obj.Trkland.cingulum.out.clineFA_lh_HDorff,'file') == 0
                    display('Executing centerline_lh for cingulum_hDorff...');
                    obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_centerline(rotrk_read(obj.Trkland.cingulum.out.clean_trks_lh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'cingulum_lh_cleantrimmed'), 'hausdorff');
                    %Adding scalars:
                    obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.cingulum_clineHDorff_lh.header,obj.Trkland.Trks.cingulum_clineHDorff_lh.sstr,obj.Trkland.cingulum.out.clineFA_lh_HDorff);
                end
                
                %GET THE DATA VALUES NOW:
                %Unclean TRKS:
                if  exist(obj.Trkland.cingulum.out.trk_lh,'file') ~= 0
                    obj.Trkland.cingulum.data.lh_unclean_vol = obj.Trkland.Trks.raw_cingulum_lh.num_uvox;
                    obj.Trkland.cingulum.data.lh_unclean_FA = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.lh_unclean_RD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.lh_unclean_AxD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.lh_unclean_MD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,7));
                else
                    
                    obj.Trkland.cingulum.data.lh_unclean_vol = [];
                    obj.Trkland.cingulum.data.lh_unclean_FA = [];
                    obj.Trkland.cingulum.data.lh_unclean_RD = [];
                    obj.Trkland.cingulum.data.lh_unclean_AxD = [];
                    obj.Trkland.cingulum.data.lh_unclean_MD = [];
                    
                    %Fill dependency trks to []:
                    obj.remove_trkland_fields('raw_cingulum_lh')
                    obj.remove_trkland_fields('cingulum_trimmed_lh')
                    obj.remove_trkland_fields('cingulum_cleantrimmed_lh')
                    obj.remove_trkland_fields('cingulum_clinehighFA_lh')
                    obj.remove_trkland_fields('cingulum_clineHDorff_lh')
                end
                %Trimmed_clean_cingulum:
                if numel(obj.Trkland.Trks.cingulum_trimmed_lh.sstr) ~= 0
                    obj.Trkland.cingulum.data.lh_trimmedclean_vol = obj.Trkland.Trks.cingulum_trimmed_lh.num_uvox;
                    obj.Trkland.cingulum.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.lh_trimmedclean_vol = [];
                    obj.Trkland.cingulum.data.lh_trimmedclean_FA = [];
                    obj.Trkland.cingulum.data.lh_trimmedclean_RD = [];
                    obj.Trkland.cingulum.data.lh_trimmedclean_AxD = [];
                    obj.Trkland.cingulum.data.lh_trimmedclean_MD = [];
                end
                
                %Clean_cingulum:
                if numel(obj.Trkland.Trks.cingulum_cleantrimmed_lh.sstr) ~= 0
                    obj.Trkland.cingulum.data.lh_clean_vol = obj.Trkland.Trks.cingulum_cleantrimmed_lh.num_uvox;
                    obj.Trkland.cingulum.data.lh_clean_FA = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.lh_clean_RD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.lh_clean_AxD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.lh_clean_MD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.lh_clean_vol = [];
                    obj.Trkland.cingulum.data.lh_clean_FA = [];
                    obj.Trkland.cingulum.data.lh_clean_RD = [];
                    obj.Trkland.cingulum.data.lh_clean_AxD = [];
                    obj.Trkland.cingulum.data.lh_clean_MD = [];
                end
                
                %Cline_HighFA
                if numel(obj.Trkland.Trks.cingulum_clinehighFA_lh.sstr) ~= 0
                    obj.Trkland.cingulum.data.lh_cline_HighFA_vol = obj.Trkland.Trks.cingulum_clinehighFA_lh.num_uvox;
                    obj.Trkland.cingulum.data.lh_cline_length_highFA = obj.Trkland.Trks.cingulum_clinehighFA_lh.maxsstrlen;
                    obj.Trkland.cingulum.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.lh_cline_HighFA_vol = [];
                    obj.Trkland.cingulum.data.lh_cline_length_highFA= [];
                    obj.Trkland.cingulum.data.lh_cline_FA_highFA = [];
                    obj.Trkland.cingulum.data.lh_cline_RD_highFA = [];
                    obj.Trkland.cingulum.data.lh_cline_AxD_highFA = [];
                    obj.Trkland.cingulum.data.lh_cline_MD_highFA = [];
                end
                
                %Cline_HDorff
                if numel(obj.Trkland.Trks.cingulum_clineHDorff_lh.sstr) ~= 0
                    obj.Trkland.cingulum.data.lh_cline_HDorff_vol = obj.Trkland.Trks.cingulum_clineHDorff_lh.num_uvox;
                    obj.Trkland.cingulum.data.lh_cline_length_HDorff = obj.Trkland.Trks.cingulum_clineHDorff_lh.maxsstrlen;
                    obj.Trkland.cingulum.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.lh_cline_HDorff_vol = [];
                    obj.Trkland.cingulum.data.lh_cline_length_HDorff= [];
                    obj.Trkland.cingulum.data.lh_cline_FA_HDorff = [];
                    obj.Trkland.cingulum.data.lh_cline_RD_HDorff = [];
                    obj.Trkland.cingulum.data.lh_cline_AxD_HDorff = [];
                    obj.Trkland.cingulum.data.lh_cline_MD_HDorff = [];
                end
            end
            
            %RIGHT SIDE:
            for tohide=1:1
                if exist(obj.Trkland.cingulum.out.clean_trks_rh ,'file') == 0 && exist(obj.Trkland.cingulum.out.trk_rh,'file') ~= 0
                    
                    clear obj.Trkland.Trks.raw_cingulum_rh obj.Trkland.Trks.cingulum_cleantrimmed_rh obj.Trkland.Trks.cingulum_clineinit_rh
                    obj.Trkland.Trks.raw_cingulum_rh = rotrk_read(obj.Trkland.cingulum.out.trk_rh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'cingulum_rh_cleantrimmed');
                    %add Scalars
                    obj.Trkland.Trks.raw_cingulum_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.raw_cingulum_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.raw_cingulum_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.raw_cingulum_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_cingulum_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    
                    %Trim tracts here (checking if trimming occurs, mainly for manual edits):
                    if exist(obj.Trkland.cingulum.out.clean_trkstrimmed_rh,'file') == 0
                        obj.Trkland.Trks.cingulum_trimmed_rh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_cingulum_rh, ...
                            [ {obj.Trkland.cingulum.in.seed_postcing_rh}  {obj.Trkland.cingulum.in.rostantcing_rh}  ], 'cingulum_rh');
                    else
                        obj.Trkland.Trks.cingulum_trimmed_rh = rotrk_read(obj.Trkland.cingulum.out.clean_trkstrimmed_rh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'cingulum_rh_cleantrimmed' );
                        obj.Trkland.Trks.cingulum_trimmed_rh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_rh  ,obj.Params.Dtifit.out.FA{end} , 'FA');
                        obj.Trkland.Trks.cingulum_trimmed_rh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                        obj.Trkland.Trks.cingulum_trimmed_rh = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                        obj.Trkland.Trks.cingulum_trimmed_rh  = rotrk_add_sc(obj.Trkland.Trks.cingulum_trimmed_rh  ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    end
                    %Select the HDorff centerline(first pass)
                    obj.Trkland.Trks.cingulum_clineinit_rh= rotrk_centerline(obj.Trkland.Trks.cingulum_trimmed_rh,'hausdorff');
                    %Clean up based on normality of hausdorff distance
                    obj.Trkland.Trks.cingulum_cleantrimmed_rh = rotrk_rm_byHDorff(obj.Trkland.Trks.cingulum_clineinit_rh, obj.Trkland.Trks.cingulum_trimmed_rh,obj.Trkland.Trks.cingulum_trimmed_rh);
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.cingulum_cleantrimmed_rh.header,obj.Trkland.Trks.cingulum_cleantrimmed_rh.sstr,obj.Trkland.cingulum.out.clean_trks_rh )
                    rotrk_write(obj.Trkland.Trks.cingulum_trimmed_rh.header,obj.Trkland.Trks.cingulum_trimmed_rh.sstr,obj.Trkland.cingulum.out.clean_trkstrimmed_rh);
                    wasRun=true;
                    obj.UpdateHist(obj.Trkland.cingulum,'trkland_cingulum', obj.Trkland.cingulum.out.clineFA_rh_highFA,wasRun);
                end
                
                %Now that the TRK is clean, lets get the high_FA and get the centerline:
                %Pick centerline based on high_sc and FA:
                %High_FA_centerline:
                if exist(obj.Trkland.cingulum.out.clineFA_rh_highFA,'file')==0
                    display('Executing centerline_rh for cingulum_highFA...');
                    temp_clean_trk_rh=rotrk_read(obj.Trkland.cingulum.out.clean_trks_rh ,obj.sessionname,obj.Params.Dtifit.out.FA{end},'cingulum_rh_cleantrimmed');
                    temp_clean_trk_rh=rotrk_add_sc(temp_clean_trk_rh,obj.Params.Dtifit.out.FA{end},'FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_rh= rotrk_centerline(temp_clean_trk_rh, 'high_sc','FA');
                    %Adding scalars:
                    obj.Trkland.Trks.cingulum_clinehighFA_rh = rotrk_centerline(temp_clean_trk_rh, 'high_sc','FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.cingulum_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.cingulum_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.cingulum_clinehighFA_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clinehighFA_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.cingulum_clinehighFA_rh.header,obj.Trkland.Trks.cingulum_clinehighFA_rh.sstr,obj.Trkland.cingulum.out.clineFA_rh_highFA );
                end
                %HDorff_centerline:
                if exist(obj.Trkland.cingulum.out.clineFA_rh_HDorff,'file')==0
                    display('Executing centerline_rh for cingulum_hDorff...');
                    obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_centerline(rotrk_read(obj.Trkland.cingulum.out.clean_trks_rh,obj.sessionname,obj.Params.Dtifit.out.FA{end},'cingulum_rh_cleantrimmed'), 'hausdorff');
                    %Adding scalars:
                    obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
                    obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
                    obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
                    obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_add_sc(  obj.Trkland.Trks.cingulum_clineHDorff_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
                    %save trks:
                    rotrk_write(obj.Trkland.Trks.cingulum_clineHDorff_rh.header,obj.Trkland.Trks.cingulum_clineHDorff_rh.sstr,obj.Trkland.cingulum.out.clineFA_rh_HDorff);
                end
                
                %GET THE DATA VALUES NOW:
                if exist(obj.Trkland.cingulum.out.trk_rh,'file') ~= 0
                    %Unclean TRKS:
                    obj.Trkland.cingulum.data.rh_unclean_vol = obj.Trkland.Trks.raw_cingulum_rh.num_uvox;
                    obj.Trkland.cingulum.data.rh_unclean_FA = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.rh_unclean_RD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.rh_unclean_AxD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.rh_unclean_MD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.rh_unclean_vol =  [];
                    obj.Trkland.cingulum.data.rh_unclean_FA = [];
                    obj.Trkland.cingulum.data.rh_unclean_RD = [];
                    obj.Trkland.cingulum.data.rh_unclean_AxD = [];
                    obj.Trkland.cingulum.data.rh_unclean_MD = [];
                    
                    %Fill dependency trks to []:
                    obj.remove_trkland_fields('raw_cingulum_rh')
                    obj.remove_trkland_fields('cingulum_trimmed_rh')
                    obj.remove_trkland_fields('cingulum_cleantrimmed_rh')
                    obj.remove_trkland_fields('cingulum_clinehighFA_rh')
                    obj.remove_trkland_fields('cingulum_clineHDorff_rh')
                end
                
                %Trimmed_clean_cingulum:
                if numel(obj.Trkland.Trks.cingulum_trimmed_rh.sstr) ~= 0
                    obj.Trkland.cingulum.data.rh_trimmedclean_vol = obj.Trkland.Trks.cingulum_trimmed_rh.num_uvox;
                    obj.Trkland.cingulum.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.rh_trimmedclean_vol = [];
                    obj.Trkland.cingulum.data.rh_trimmedclean_FA = [];
                    obj.Trkland.cingulum.data.rh_trimmedclean_RD = [];
                    obj.Trkland.cingulum.data.rh_trimmedclean_AxD = [];
                    obj.Trkland.cingulum.data.rh_trimmedclean_MD = [];
                end
                
                %Clean_cingulum:
                if numel(obj.Trkland.Trks.cingulum_cleantrimmed_rh.sstr) ~= 0
                    obj.Trkland.cingulum.data.rh_clean_vol = obj.Trkland.Trks.cingulum_cleantrimmed_rh.num_uvox;
                    obj.Trkland.cingulum.data.rh_clean_FA = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.rh_clean_RD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.rh_clean_AxD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.rh_clean_MD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.rh_clean_vol = [];
                    obj.Trkland.cingulum.data.rh_clean_FA = [];
                    obj.Trkland.cingulum.data.rh_clean_RD = [];
                    obj.Trkland.cingulum.data.rh_clean_AxD = [];
                    obj.Trkland.cingulum.data.rh_clean_MD = [];
                end
                
                %Cline_HighFA
                if numel(obj.Trkland.Trks.cingulum_clinehighFA_rh.sstr) ~= 0
                    obj.Trkland.cingulum.data.rh_cline_HighFA_vol = obj.Trkland.Trks.cingulum_clinehighFA_rh.num_uvox;
                    obj.Trkland.cingulum.data.rh_cline_length_highFA=obj.Trkland.Trks.cingulum_clinehighFA_rh.maxsstrlen;
                    obj.Trkland.cingulum.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.rh_cline_HighFA_vol = [];
                    obj.Trkland.cingulum.data.rh_cline_length_highFA= [];
                    obj.Trkland.cingulum.data.rh_cline_FA_highFA = [];
                    obj.Trkland.cingulum.data.rh_cline_RD_highFA = [];
                    obj.Trkland.cingulum.data.rh_cline_AxD_highFA = [];
                    obj.Trkland.cingulum.data.rh_cline_MD_highFA = [];
                end
                
                %Cline_HDorff
                if numel(obj.Trkland.Trks.cingulum_clineHDorff_rh.sstr) ~= 0
                    obj.Trkland.cingulum.data.rh_cline_HDorff_vol = obj.Trkland.Trks.cingulum_clineHDorff_rh.num_uvox;
                    obj.Trkland.cingulum.data.rh_cline_length_HDorff=obj.Trkland.Trks.cingulum_clineHDorff_rh.maxsstrlen;
                    obj.Trkland.cingulum.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,4));
                    obj.Trkland.cingulum.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,5));
                    obj.Trkland.cingulum.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,6));
                    obj.Trkland.cingulum.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,7));
                else
                    obj.Trkland.cingulum.data.rh_cline_HDorff_vol = [];
                    obj.Trkland.cingulum.data.rh_cline_length_HDorff= [];
                    obj.Trkland.cingulum.data.rh_cline_FA_HDorff = [];
                    obj.Trkland.cingulum.data.rh_cline_RD_HDorff = [];
                    obj.Trkland.cingulum.data.rh_cline_AxD_HDorff = [];
                    obj.Trkland.cingulum.data.rh_cline_MD_HDorff = [];
                end
            end
        end
        
        function obj = trkland_atr(obj)
            display( 'trkland_atr implementation has not been. Please check before using it (improvements were given to trimmmed_clean and ordeing so check. Dat stamped 10162017!!!');
            display('Skipping trkland_atr() ...')
            
%             wasRun = false;
%             %Create trkland directory (if doesn't exist)
%             exec_cmd = [ 'mkdir -p ' obj.Trkland.root ];
%             obj.RunBash(exec_cmd);
%             outpath=obj.Trkland.root;
%             
%             %ROIs/SEEDs PREPARATION
%             
%             obj.Trkland.atr.in.roi_antroscing_lh = [ obj.Trkland.root 'atr_roi_antrostralcingulateDil1_lh.nii.gz' ] ;
%             obj.Trkland.atr.in.roi_antroscing_rh = [ obj.Trkland.root 'atr_roi_antrostralcingulateDil1_rh.nii.gz' ] ;
%             
%             obj.Trkland.atr.in.seed_thalamus_lh = [ obj.Trkland.root 'atr_seed_thalamusDil1_lh.nii.gz' ] ;
%             obj.Trkland.atr.in.seed_thalamus_rh = [ obj.Trkland.root 'atr_seed_thalamusDil1_rh.nii.gz' ] ;
%             for tohide=1:1
%                 %anterior cingulate dilation:
%                 if exist(obj.Trkland.atr.in.roi_antroscing_lh , 'file') == 0
%                     fprintf('\ntrkland_atr(): Working on anterior cingulate_lh...')
%                     exec_cmd = ['fslmaths '  obj.Trkland.atr.in.rostantcing_lh ...
%                         ' -dilM ' obj.Trkland.atr.in.roi_antroscing_lh  ];
%                     obj.RunBash(exec_cmd);
%                     fprintf('...done \n');
%                 end
%                 
%                 if exist(obj.Trkland.atr.in.roi_antroscing_rh , 'file') == 0
%                     fprintf('\ntrkland_atr(): Working on anterior cingulate_rh...')
%                     exec_cmd = ['fslmaths '  obj.Trkland.atr.in.rostantcing_rh ...
%                         ' -dilM ' obj.Trkland.atr.in.roi_antroscing_rh  ];
%                     obj.RunBash(exec_cmd);
%                     fprintf('...done \n');
%                 end
%                 
%                 %Thalamus dilation
%                 if exist(obj.Trkland.atr.in.seed_thalamus_lh , 'file') == 0
%                     fprintf('\n Working on posterior cingulate_lh...')
%                     exec_cmd = ['fslmaths ' obj.Trkland.atr.in.thalamus_lh ...
%                         ' -dilM ' obj.Trkland.atr.in.seed_thalamus_lh  ];
%                     obj.RunBash(exec_cmd);
%                     fprintf('...done \n');
%                 end
%                 if exist(obj.Trkland.atr.in.seed_thalamus_rh , 'file') == 0
%                     fprintf('\n Working on posterior cingulate_rh...')
%                     exec_cmd = ['fslmaths ' obj.Trkland.atr.in.thalamus_rh ...
%                         ' -dilM ' obj.Trkland.atr.in.seed_thalamus_rh  ];
%                     obj.RunBash(exec_cmd);
%                     fprintf('...done \n');
%                 end
%                 
%             end
%             
%             %CLEANUP OF THE TRACTS
%             for tohide=1:1
%                 obj.Trkland.atr.out.clean_trkstrimmed_lh = [ obj.Trkland.root  'trkk_atr_trimmed_lh.trk.gz'];
%                 obj.Trkland.atr.out.clean_trkstrimmed_rh = [ obj.Trkland.root  'trkk_atr_trimmed_rh.trk.gz'];
%                 
%                 obj.Trkland.atr.out.clean_trks_lh = [ obj.Trkland.root  'trkk_atr_trimmedclean_lh.trk.gz'];
%                 obj.Trkland.atr.out.clean_trks_rh = [ obj.Trkland.root  'trkk_atr_trimmedclean_rh.trk.gz'];
%                 
%                 obj.Trkland.atr.out.clineFA_lh_highFA = [ obj.Trkland.root  'cline_atr_highFAlh.trk.gz'];
%                 obj.Trkland.atr.out.clineFA_lh_HDorff = [ obj.Trkland.root  'cline_atr_HDorfflh.trk.gz'];
%                 obj.Trkland.atr.out.clineFA_rh_highFA = [ obj.Trkland.root  'cline_atr_highFArh.trk.gz'];
%                 obj.Trkland.atr.out.clineFA_rh_HDorff = [ obj.Trkland.root  'cline_atr_HDorffrh.trk.gz'];
%                 
%                 obj.Trkland.atr.cingulum.QCfile_lh = [outpath 'QC_atr_lh.flag'] ;
%                 obj.Trkland.atr.cingulum.QCfile_rh = [outpath 'QC_atr_rh.flag'] ;
%                 obj.Trkland.atr.cingulum.QCfile_bil = [outpath 'QC_atr_bil.flag'] ;
%                 
%             end
%             
%             %TRACKING STARS HERE:
%             for tohide=1:1
%                 obj.Trkland.atr.out.trk_lh = [ obj.Trkland.root  'trkk_atr_lh.trk.gz'];
%                 obj.Trkland.atr.out.trk_rh = [ obj.Trkland.root  'trkk_atr_rh.trk.gz'];
%                 if exist(obj.Trkland.atr.cingulum.QCfile_bil,'file') == 0
%                     %Left side trking:
%                     if exist( obj.Trkland.atr.cingulum.QCfile_lh,'file')==0
%                         if exist(obj.Trkland.atr.out.trk_lh,'file') == 0
%                             exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
%                                 ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
%                                 ' --seed=' obj.Trkland.atr.in.seed_thalamus_lh ' --roi=' obj.Trkland.atr.in.roi_antroscing_lh ...
%                                 ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
%                                 ' --output=' obj.Trkland.atr.out.trk_lh ];
%                             for dd=1:4 %trying 4 times to get a trk. If not, quit!
%                                 if exist(obj.Trkland.atr.out.trk_lh,'file') == 0
%                                     obj.RunBash(exec_cmd,144);
%                                 end
%                             end
%                             wasRun=true;
%                             obj.UpdateHist(obj.Trkland.atr,'trkland_atr', obj.Trkland.atr.out.trk_lh,wasRun);
%                         end
%                     else
%                         display('QC_flag_rh found in trkland_atr. Skipping and removing data points...')
%                         RefreshFields(obj,'atr','rh')
%                     end
%                     %Right side trking:
%                     if exist( obj.Trkland.atr.cingulum.QCfile_rh,'file')==0
%                         if exist(obj.Trkland.atr.out.trk_rh,'file') == 0
%                             exec_cmd = ['dsi_studio_run --action=trk --source=' obj.Trkland.fx.in.fib ...
%                                 ' --seed_count=20000 --smoothing=0.01 --method=0 --interpolation=0 --thread_count=10' ...
%                                 ' --seed=' obj.Trkland.atr.in.seed_thalamus_rh ' --roi=' obj.Trkland.atr.in.roi_antroscing_rh ...
%                                 ' --step_size=1 --turning_angle=40 --min_length=110 --max_length=250 ' ...
%                                 ' --output=' obj.Trkland.atr.out.trk_rh ];
%                             
%                             for dd=1:4 %trying 4 times to get a trk. If not, quit!
%                                 if exist(obj.Trkland.atr.out.trk_rh,'file') == 0
%                                     obj.RunBash(exec_cmd,144);
%                                 end
%                             end
%                             wasRun=true;
%                             obj.UpdateHist(obj.Trkland.atr,'trkland_atr', obj.Trkland.atr.out.trk_rh,wasRun);
%                             
%                         end
%                     else
%                         display('QC_flag_lh found in trkland_atr. Skipping and removing data points...')
%                         RefreshFields(obj,'atr','rh')
%                     end
%                 else
%                     display('QC_flag_bil found in trkland_atr. Skipping and removing data points...')
%                     RefreshFields(obj,'atr','bil')
%                 end
%             end
%             
%             
%             %LEFT SIDE:
%             for tohide=1:1
%                 if exist(obj.Trkland.atr.out.clean_trks_lh ,'file') == 0 && exist(obj.Trkland.atr.out.trk_lh,'file') ~= 0
%                     try
%                         obj.Trkland.Trks.raw_atr_lh = rotrk_read(obj.Trkland.atr.out.trk_lh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'atr_lh');
%                         %add Scalars:
%                         obj.Trkland.Trks.raw_atr_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_lh ,obj.Params.Dtifit.out.FA{end} , 'FA');
%                         obj.Trkland.Trks.raw_atr_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
%                         obj.Trkland.Trks.raw_atr_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
%                         obj.Trkland.Trks.raw_atr_lh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_lh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
%                         %Trim tracts here:
%                         obj.Trkland.Trks.atr_trimmed_lh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_atr_lh, ...
%                             [ { obj.Trkland.atr.in.thalamus_lh}  {obj.Trkland.atr.in.rostantcing_lh}  ], 'atr_lh');
%                         
%                         %Select the HDorff centerline(first pass):
%                         obj.Trkland.Trks.atr_clineinit_lh= rotrk_centerline(obj.Trkland.Trks.atr_trimmed_lh,'hausdorff');
%                         %Clean up based on normality of hausdorff distance:
%                         obj.Trkland.Trks.atr_cleantrimmed_lh = rotrk_rm_byHDorff(obj.Trkland.Trks.atr_clineinit_lh, obj.Trkland.Trks.atr_trimmed_lh,obj.Trkland.Trks.atr_trimmed_lh);
%                         %%%obj.Trkland.Trks.atr_cleantrimmed_lh = rotrk_rm_bylen(temp_fx_lh_cline, temp_fx_lh,temp_fx_lh);
%                         %Now that the TRK is clean, lets get the high_FA and get the centerline:
%                         %Pick centerline based on high_sc and FA:
%                         obj.Trkland.Trks.atr_clinehighFA_lh= rotrk_centerline(obj.Trkland.Trks.atr_cleantrimmed_lh, 'high_sc','FA');
%                         obj.Trkland.Trks.atr_clineHDorff_lh = rotrk_centerline(obj.Trkland.Trks.atr_cleantrimmed_lh, 'hausdorff');
%                         %save trks:
%                         rotrk_write(obj.Trkland.Trks.atr_trimmed_lh.header,obj.Trkland.Trks.atr_trimmed_lh.sstr,obj.Trkland.atr.out.clean_trkstrimmed_lh);
%                         rotrk_write(obj.Trkland.Trks.atr_cleantrimmed_lh.header,obj.Trkland.Trks.atr_cleantrimmed_lh.sstr,obj.Trkland.atr.out.clean_trks_lh );
%                         rotrk_write(obj.Trkland.Trks.atr_clinehighFA_lh.header,obj.Trkland.Trks.atr_clinehighFA_lh.sstr,obj.Trkland.atr.out.clineFA_lh_highFA );
%                         rotrk_write(obj.Trkland.Trks.atr_clineHDorff_lh.header,obj.Trkland.Trks.atr_clineHDorff_lh.sstr,obj.Trkland.atr.out.clineFA_lh_HDorff);
%                         wasRun=true;
%                         obj.UpdateHist(obj.Trkland.atr,'trkland_atr', obj.Trkland.atr.out.clineFA_lh_highFA,wasRun);
%                         %Get volume data of unclean/cleaned tracts:
%                         %Volume:
%                         obj.Trkland.atr.data.lh_unclean_vol = obj.Trkland.Trks.raw_atr_lh.num_uvox;
%                         obj.Trkland.atr.data.lh_trimmedclean_vol = obj.Trkland.Trks.atr_trimmed_lh.num_uvox;
%                         obj.Trkland.atr.data.lh_clean_vol = obj.Trkland.Trks.atr_cleantrimmed_lh.num_uvox;
%                         obj.Trkland.atr.data.lh_cline_HighFA_vol = obj.Trkland.Trks.atr_clinehighFA_lh.num_uvox;
%                         obj.Trkland.atr.data.lh_cline_HDorff_vol = obj.Trkland.Trks.atr_clineHDorff_lh.num_uvox;
%                         %METRICS DATA NOW
%                         %unclean_atr_lh
%                         obj.Trkland.atr.data.lh_unclean_FA = mean(obj.Trkland.Trks.raw_atr_lh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.lh_unclean_RD = mean(obj.Trkland.Trks.raw_atr_lh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.lh_unclean_AxD = mean(obj.Trkland.Trks.raw_atr_lh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.lh_unclean_MD = mean(obj.Trkland.Trks.raw_atr_lh.unique_voxels(:,7));
%                         %Clean_atr:
%                         obj.Trkland.atr.data.lh_clean_FA = mean(obj.Trkland.Trks.atr_cleantrimmed_lh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.lh_clean_RD = mean(obj.Trkland.Trks.atr_cleantrimmed_lh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.lh_clean_AxD = mean(obj.Trkland.Trks.atr_cleantrimmed_lh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.lh_clean_MD = mean(obj.Trkland.Trks.atr_cleantrimmed_lh.unique_voxels(:,7));
%                         %Trimmed_clean_atr:
%                         obj.Trkland.atr.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.atr_trimmed_lh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.atr_trimmed_lh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.atr_trimmed_lh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.atr_trimmed_lh.unique_voxels(:,7));
%                         %Cline_HighFA
%                         obj.Trkland.atr.data.lh_cline_length_highFA=obj.Trkland.Trks.atr_clinehighFA_lh.maxsstrlen;
%                         obj.Trkland.atr.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_lh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_lh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_lh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_lh.unique_voxels(:,7));
%                         %Cline_HDorff
%                         obj.Trkland.atr.data.lh_cline_length_HDorff=obj.Trkland.Trks.atr_clineHDorff_lh.maxsstrlen;
%                         obj.Trkland.atr.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_lh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_lh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_lh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_lh.unique_voxels(:,7));
%                     catch
%                         warning('No cleanup for atr_lh was finished correctly. Maybe the specific tract cant be reached. Check raw tracts...')
%                     end
%                 end
%             end
%             
%             %RIGHT SIDE:
%             for tohide=1:1
%                 if exist(obj.Trkland.atr.out.clean_trks_rh ,'file') == 0 && exist(obj.Trkland.atr.out.trk_rh,'file') ~= 0
%                     try
%                         clear obj.Trkland.Trks.raw_atr_rh obj.Trkland.Trks.atr_cleantrimmed_rh obj.Trkland.Trks.atr_clineinit_rh
%                         obj.Trkland.Trks.raw_atr_rh = rotrk_read(obj.Trkland.atr.out.trk_rh, obj.sessionname, obj.Params.Dtifit.out.FA{end}, 'atr_rh');
%                         %add Scalars
%                         obj.Trkland.Trks.raw_atr_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_rh ,obj.Params.Dtifit.out.FA{end} , 'FA');
%                         obj.Trkland.Trks.raw_atr_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','RD') , 'RD');
%                         obj.Trkland.Trks.raw_atr_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','AxD') , 'AxD');
%                         obj.Trkland.Trks.raw_atr_rh = rotrk_add_sc(  obj.Trkland.Trks.raw_atr_rh ,strrep(obj.Params.Dtifit.out.FA{end},'FA','MD') , 'MD');
%                         %Trim tracts here:
%                         obj.Trkland.Trks.atr_trimmed_rh = rotrk_trimmedbyTOI(obj.Trkland.Trks.raw_atr_rh, ...
%                             [  { obj.Trkland.atr.in.thalamus_rh}  {obj.Trkland.atr.in.rostantcing_rh}    ], 'atr_rh');
%                         %
%                         %Select the HDorff centerline(first pass)
%                         obj.Trkland.Trks.atr_clineinit_rh= rotrk_centerline(obj.Trkland.Trks.atr_trimmed_rh,'hausdorff');
%                         %Clean up based on normality of hausdorff distance
%                         obj.Trkland.Trks.atr_cleantrimmed_rh = rotrk_rm_byHDorff(obj.Trkland.Trks.atr_clineinit_rh, obj.Trkland.Trks.atr_trimmed_rh,obj.Trkland.Trks.atr_trimmed_rh);
%                         %%%%obj.Trkland.Trks.atr_cleantrimmed_rh = rotrk_rm_bylen(temp_fx_rh_cline, temp_fx_rh,temp_fx_rh);
%                         %Now that the TRK is clean, lets get the high_FA and get the centerline:
%                         %Pick centerline based on high_sc and FA:
%                         obj.Trkland.Trks.atr_clinehighFA_rh= rotrk_centerline(obj.Trkland.Trks.atr_cleantrimmed_rh, 'high_sc','FA');
%                         obj.Trkland.Trks.atr_clineHDorff_rh = rotrk_centerline(obj.Trkland.Trks.atr_cleantrimmed_rh, 'hausdorff');
%                         %save trks:
%                         rotrk_write(obj.Trkland.Trks.atr_cleantrimmed_rh.header,obj.Trkland.Trks.atr_cleantrimmed_rh.sstr,obj.Trkland.atr.out.clean_trks_rh )
%                         rotrk_write(obj.Trkland.Trks.atr_trimmed_rh.header,obj.Trkland.Trks.atr_trimmed_rh.sstr,obj.Trkland.atr.out.clean_trkstrimmed_rh);
%                         rotrk_write(obj.Trkland.Trks.atr_clinehighFA_rh.header,obj.Trkland.Trks.atr_clinehighFA_rh.sstr,obj.Trkland.atr.out.clineFA_rh_highFA )
%                         rotrk_write(obj.Trkland.Trks.atr_clineHDorff_rh.header,obj.Trkland.Trks.atr_clineHDorff_rh.sstr,obj.Trkland.atr.out.clineFA_rh_HDorff)
%                         wasRun=true;
%                         obj.UpdateHist(obj.Trkland.atr,'trkland_atr', obj.Trkland.atr.out.clineFA_rh_highFA,wasRun);
%                         %Get volume data of unclean/cleaned tracts:
%                         %Volume:
%                         obj.Trkland.atr.data.rh_unclean_vol = obj.Trkland.Trks.raw_atr_rh.num_uvox;
%                         obj.Trkland.atr.data.rh_clean_vol = obj.Trkland.Trks.atr_cleantrimmed_rh.num_uvox;
%                         obj.Trkland.atr.data.rh_trimmedclean_vol = obj.Trkland.Trks.atr_trimmed_rh.num_uvox;
%                         obj.Trkland.atr.data.rh_cline_HighFA_vol = obj.Trkland.Trks.atr_clinehighFA_rh.num_uvox;
%                         obj.Trkland.atr.data.rh_cline_HDorff_vol = obj.Trkland.Trks.atr_clineHDorff_rh.num_uvox;
%                         %METRICS DATA NOW
%                         %unclean_atr_rh
%                         obj.Trkland.atr.data.rh_unclean_FA = mean(obj.Trkland.Trks.raw_atr_rh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.rh_unclean_RD = mean(obj.Trkland.Trks.raw_atr_rh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.rh_unclean_AxD = mean(obj.Trkland.Trks.raw_atr_rh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.rh_unclean_MD = mean(obj.Trkland.Trks.raw_atr_rh.unique_voxels(:,7));
%                         %Clean_atr:
%                         obj.Trkland.atr.data.rh_clean_FA = mean(obj.Trkland.Trks.atr_cleantrimmed_rh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.rh_clean_RD = mean(obj.Trkland.Trks.atr_cleantrimmed_rh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.rh_clean_AxD = mean(obj.Trkland.Trks.atr_cleantrimmed_rh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.rh_clean_MD = mean(obj.Trkland.Trks.atr_cleantrimmed_rh.unique_voxels(:,7));
%                         %Trimmed_clean_atr:
%                         obj.Trkland.atr.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.atr_trimmed_rh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.atr_trimmed_rh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.atr_trimmed_rh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.atr_trimmed_rh.unique_voxels(:,7));
%                         %Cline_HighFA
%                         obj.Trkland.atr.data.rh_cline_length_highFA=obj.Trkland.Trks.atr_clinehighFA_rh.maxsstrlen;
%                         obj.Trkland.atr.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_rh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_rh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_rh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.atr_clinehighFA_rh.unique_voxels(:,7));
%                         %Cline_HDorff
%                         obj.Trkland.atr.data.rh_cline_length_HDorff=obj.Trkland.Trks.atr_clineHDorff_rh.maxsstrlen;
%                         obj.Trkland.atr.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_rh.unique_voxels(:,4));
%                         obj.Trkland.atr.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_rh.unique_voxels(:,5));
%                         obj.Trkland.atr.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_rh.unique_voxels(:,6));
%                         obj.Trkland.atr.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.atr_clineHDorff_rh.unique_voxels(:,7));
%                     catch
%                         warning('No cleanup for atr_rh was finished correctly. Maybe the specific tract cant be reached? Check raw tracts...')
%                     end
%                 end
%             end
        end
        
        function obj = proc_tracx2thal11(obj)
            % Make sure you
            wasRun=false;
            %Creating Proc specific Directory:
            obj.Params.tracx_thal2ctx11.in_dir=obj.getPath(obj.Params.tracx_thal2ctx11.in.bedp_dir,obj.Params.tracx_thal2ctx11.in.movefiles);
            
            
            %PREPARING 10 CTX SEGMENTATIONS (per L/H hemispheres)
            %check if the dependency list exists...
            if exist(obj.Params.tracx_thal2ctx11.in.prep_segs_list,'file') == 0
                error(['\n proc_tracx2thal11(): ' obj.Params.tracx_thal2ctx11.in.prep_segs_list ' does not exit' ]);
            else
                tmp_readTXT=fileread([ obj.Params.tracx_thal2ctx11.in.prep_segs_list ]);
                
            end
            obj.Params.tracx_thal2ctx11.in.temp_list=textscan(tmp_readTXT,' %s %s %s %s %s %s %s %s %s %s %s');
            
            %Creating segs dir
            obj.Params.tracx_thal2ctx11.in.segs_dir = [ obj.Params.tracx_thal2ctx11.in_dir 'segs' filesep ]
            exec_cmd = [ 'mkdir -p ' obj.Params.tracx_thal2ctx11.in.segs_dir ];
            obj.RunBash(exec_cmd);
            
            
            %Creating segmentations...
            obj.Params.tracx_thal2ctx11.in.seg_list = '' ;
            for ii=1:size(obj.Params.tracx_thal2ctx11.in.temp_list,2)
                %For loop on every column that makes each iteration
                temp_seg_L = ' '; temp_seg_R = ' ';
                cur_FSseg_L = ' ' ; cur_FSseg_R =  ' ' ;
                for jj=2:size(obj.Params.tracx_thal2ctx11.in.temp_list{ii},1)
                    %check if current segmentation exists ("NA" no character), if not quit and
                    %throw error!
                    if ~strcmp(obj.Params.tracx_thal2ctx11.in.temp_list{ii}{jj},'NA')
                        cur_FSseg_L = [ obj.Params.tracx_thal2ctx11.in.FSaparc_dir ...
                            'dwi_ctx-lh-'  obj.Params.tracx_thal2ctx11.in.temp_list{ii}{jj} '.nii.gz'];
                        %check lh here:
                        if exist(cur_FSseg_L,'file') == 0
                            error(['In proc_tracx2thal11(): ' cur_FSseg_L ' does not exist. Needed for merge segmentations']);
                        end
                        
                        cur_FSseg_R = [ obj.Params.tracx_thal2ctx11.in.FSaparc_dir ...
                            'dwi_ctx-rh-'  obj.Params.tracx_thal2ctx11.in.temp_list{ii}{jj} '.nii.gz' ];
                        %check rh here:
                        if exist(cur_FSseg_R,'file') == 0
                            error(['\nIn proc_tracx2thal11(): ' cur_FSseg_R ' does not exist. Needed for merge segmentations']);
                        end
                        %Initial "-add" is not needed for first element
                        if jj == 2
                            temp_seg_L = [ temp_seg_L  ' ' cur_FSseg_L ' '   ];
                            temp_seg_R = [ temp_seg_R  ' ' cur_FSseg_R ' ' ];
                        else
                            temp_seg_L = [ temp_seg_L ' -add ' cur_FSseg_L  ' ' ];
                            temp_seg_R = [ temp_seg_R ' -add ' cur_FSseg_R  ' ' ];
                        end
                    end
                end
                
                %Add all directories in a merge file using fslmaths:
                %lh:
                cur_segname_L = [ obj.Params.tracx_thal2ctx11.in.segs_dir 'lh_' ...
                    obj.Params.tracx_thal2ctx11.in.temp_list{ii}{1} '.nii.gz' ];
                exec_cmd = ['fslmaths ' temp_seg_L  cur_segname_L] ;
                %Creating a list to be saved:
                obj.Params.tracx_thal2ctx11.in.seg_list{ii} = obj.Params.tracx_thal2ctx11.in.temp_list{ii}{1};
                
                
                if exist(cur_segname_L,'file') == 0
                    display([' In proc_tracx2thal11(): merging ' cur_segname_L '...'])
                    obj.RunBash(exec_cmd);
                    display('...done');
                end
                %rh:
                cur_segname_R = [ obj.Params.tracx_thal2ctx11.in.segs_dir 'rh_' ...
                    obj.Params.tracx_thal2ctx11.in.temp_list{ii}{1} '.nii.gz' ];
                exec_cmd = ['fslmaths ' temp_seg_R cur_segname_R ] ;
                if exist(cur_segname_R) == 0
                    display([' In proc_tracx2thal11(): merging ' cur_segname_R '...'])
                    obj.RunBash(exec_cmd);
                    display('...done');
                end
            end
            
            %Creating the txt files:
            %lh:
            obj.Params.tracx_thal2ctx11.in.lh_txt =  [ obj.Params.tracx_thal2ctx11.in.segs_dir 'segs_lh.txt' ];
            exec_cmd = [ ' ls -1 ' obj.Params.tracx_thal2ctx11.in.segs_dir ...
                'lh_* > ' obj.Params.tracx_thal2ctx11.in.lh_txt  ];
            obj.RunBash(exec_cmd);
            %rh:
            obj.Params.tracx_thal2ctx11.in.rh_txt =  [ obj.Params.tracx_thal2ctx11.in.segs_dir 'segs_rh.txt' ];
            exec_cmd = [ ' ls -1 ' obj.Params.tracx_thal2ctx11.in.segs_dir ...
                'rh_* > ' obj.Params.tracx_thal2ctx11.in.rh_txt  ];
            obj.RunBash(exec_cmd);
            
            
            
            %Now ready to perform probablistic tractograhy with the
            %parameters of interest:
            %Selecingt our seeds of interes:
            obj.Params.tracx_thal2ctx11.in.lh_seed = strrep(obj.Params.tracx_thal2ctx11.in.FSaparc_dir,['_aseg' filesep'],['2009_aseg' filesep 'dwi_fs_Left-Thalamus-Proper.nii.gz']);
            obj.Params.tracx_thal2ctx11.in.rh_seed = strrep(obj.Params.tracx_thal2ctx11.in.FSaparc_dir,['_aseg' filesep'],['2009_aseg' filesep 'dwi_fs_Right-Thalamus-Proper.nii.gz']);
            
            %lh:
            exec_cmd = ['probtrackx2 -x ' obj.Params.tracx_thal2ctx11.in.lh_seed ...
                ' -l --loopcheck -c 0.2 -S 2000 --steplength=0.5 -P 5000' ...
                ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s ' ...
                ' ' obj.Params.tracx_thal2ctx11.in.bedp_dir filesep 'merged' ...
                ' -m ' obj.Params.tracx_thal2ctx11.in.bedp_dir filesep 'nodif_brain_mask' ...
                ' --dir=' obj.Params.tracx_thal2ctx11.in_dir ' --otargetpaths' ...
                ' --targetmasks=' obj.Params.tracx_thal2ctx11.in.lh_txt ' --os2t'];
            
            
            
            obj.Params.tracx_thal2ctx11.out.seed2temporal_lh = [ obj.Params.tracx_thal2ctx11.in_dir 'seeds_to_lh_temporal_3.nii.gz'] ;
            if exist(obj.Params.tracx_thal2ctx11.out.seed2temporal_lh) == 0
                obj.RunBash(exec_cmd,44);
            end
            %rh:
            exec_cmd = ['probtrackx2 -x ' obj.Params.tracx_thal2ctx11.in.rh_seed ...
                ' -l --loopcheck -c 0.2 -S 2000 --steplength=0.5 -P 5000' ...
                ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s ' ...
                ' ' obj.Params.tracx_thal2ctx11.in.bedp_dir filesep 'merged' ...
                ' -m ' obj.Params.tracx_thal2ctx11.in.bedp_dir filesep 'nodif_brain_mask' ...
                ' --dir=' obj.Params.tracx_thal2ctx11.in_dir ' --otargetpaths' ...
                ' --targetmasks=' obj.Params.tracx_thal2ctx11.in.rh_txt ' --os2t'];
            
            obj.Params.tracx_thal2ctx11.out.seed2temporal_rh = [ obj.Params.tracx_thal2ctx11.in_dir 'seeds_to_lh_temporal_3.nii.gz'] ;
            if exist(obj.Params.tracx_thal2ctx11.out.seed2temporal_rh,'file') == 0
                obj.RunBash(exec_cmd,44);
            end
            
            %Now find the biggest
            %lh:
            obj.Params.tracx_thal2ctx11.out.biggest_lh = [ obj.Params.tracx_thal2ctx11.in_dir 'lh_thal2ctx11.nii.gz'] ;
            if exist(obj.Params.tracx_thal2ctx11.out.biggest_lh,'file') == 0
                exec_cmd = ['find_the_biggest ' obj.Params.tracx_thal2ctx11.in_dir  ...
                    'seeds_to_lh_* '  obj.Params.tracx_thal2ctx11.out.biggest_lh ]
                obj.RunBash(exec_cmd,44);
            end
            %Extracting the values and uploading data
            [check_ok , tmp_lh_vals] = system(['fslstats ' ...
                obj.Params.tracx_thal2ctx11.out.biggest_lh  ' -h ' num2str([1+numel(obj.Params.tracx_thal2ctx11.in.seg_list)]) ]);
            if check_ok ~= 0
                error('In  proc_tracx2thal11(): failed to parcellate find_the_biggest for thal2ctx11_lh');
            end
            %formating the values for double characters
            floats_tmp_lh_vals=textscan(tmp_lh_vals,'%f');
            obj.Params.tracx_thal2ctx11.out.biggest_vals_lh = floats_tmp_lh_vals{1};
            
            
            %rh
            obj.Params.tracx_thal2ctx11.out.biggest_rh = [ obj.Params.tracx_thal2ctx11.in_dir 'rh_thal2ctx11.nii.gz'] ;
            if exist(obj.Params.tracx_thal2ctx11.out.biggest_rh,'file') == 0
                exec_cmd = ['find_the_biggest ' obj.Params.tracx_thal2ctx11.in_dir  ...
                    'seeds_to_rh_* '  obj.Params.tracx_thal2ctx11.out.biggest_rh ]
                obj.RunBash(exec_cmd,44);
            end
            
            
            %Extracting the values and uploading data
            [check_ok , obj.Params.tracx_thal2ctx11.out.biggest_vals_rh] = system(['fslstats ' ...
                obj.Params.tracx_thal2ctx11.out.biggest_rh  ' -h ' num2str([1+numel(obj.Params.tracx_thal2ctx11.in.seg_list)]) ])
            if check_ok ~= 0
                error('In  proc_tracx2thal11(): failed to parcellate find_the_biggest for thal2ctx11_rh' );
            end
            
            
            %Saving object:
            obj.resave
            
        end
        
        function obj = proc_tracx2papez(obj)
            % Make sure you
            wasRun=false;
            %Creating Proc specific Directory:
            obj.Params.tracx_thal2papez.in_dir=obj.getPath(obj.Params.tracx_thal2papez.in.bedp_dir,obj.Params.tracx_thal2papez.in.movefiles);
            
            
            %PREPARING 10 CTX SEGMENTATIONS (per L/H hemispheres)
            %check if the dependency list exists...
            if exist(obj.Params.tracx_thal2papez.in.prep_segs_list,'file') == 0
                error(['\n proc_tracx2papez(): ' obj.Params.tracx_thal2papez.in.prep_segs_list ' does not exit' ]);
            else
                tmp_readTXT=fileread([ obj.Params.tracx_thal2papez.in.prep_segs_list ]);
                
            end
            obj.Params.tracx_thal2papez.in.temp_list=textscan(tmp_readTXT,' %s %s %s %s %s %s %s %s ');
            
            %Creating segs dir
            obj.Params.tracx_thal2papez.in.segs_dir = [ obj.Params.tracx_thal2papez.in_dir 'segs' filesep ];
            exec_cmd = [ 'mkdir -p ' obj.Params.tracx_thal2papez.in.segs_dir ];
            obj.RunBash(exec_cmd);
            
            
            %Creating segmentations...
            obj.Params.tracx_thal2papez.in.seg_list = '' ;
            for ii=1:size(obj.Params.tracx_thal2papez.in.temp_list,2)
                %For loop on every column that makes each iteration
                temp_seg_L = ' '; temp_seg_R = ' ';
                cur_FSseg_L = ' ' ; cur_FSseg_R =  ' ' ;
                for jj=2:size(obj.Params.tracx_thal2papez.in.temp_list{ii},1)
                    %check if current segmentation exists ("NA" no character), if not quit and
                    %throw error!
                    if ~strcmp(obj.Params.tracx_thal2papez.in.temp_list{ii}{jj},'NA')
                        cur_FSseg_L = [ obj.Params.tracx_thal2papez.in.FSaparc_dir ...
                            'dwi_ctx-lh-'  obj.Params.tracx_thal2papez.in.temp_list{ii}{jj} '.nii.gz'];
                        %check lh here:
                        if exist(cur_FSseg_L,'file') == 0
                            error(['In proc_tracx2papez(): ' cur_FSseg_L ' does not exist. Needed for merge segmentations.  ij={' num2str(ii) '}{' num2str(jj) '}']);
                        end
                        
                        cur_FSseg_R = [ obj.Params.tracx_thal2papez.in.FSaparc_dir ...
                            'dwi_ctx-rh-'  obj.Params.tracx_thal2papez.in.temp_list{ii}{jj} '.nii.gz' ];
                        %check rh here:
                        if exist(cur_FSseg_R,'file') == 0
                            error(['\nIn proc_tracx2papez(): ' cur_FSseg_R ' does not exist. Needed for merge segmentations. ij={' num2str(ii) '}{' num2str(jj) '}']);
                        end
                        %Initial "-add" is not needed for first element
                        if jj == 2
                            temp_seg_L = [ temp_seg_L  ' ' cur_FSseg_L ' '   ];
                            temp_seg_R = [ temp_seg_R  ' ' cur_FSseg_R ' ' ];
                        else
                            temp_seg_L = [ temp_seg_L ' -add ' cur_FSseg_L  ' ' ];
                            temp_seg_R = [ temp_seg_R ' -add ' cur_FSseg_R  ' ' ];
                        end
                    end
                end
                
                %Add all directories in a merge file using fslmaths:
                %lh:
                cur_segname_L = [ obj.Params.tracx_thal2papez.in.segs_dir 'lh_' ...
                    obj.Params.tracx_thal2papez.in.temp_list{ii}{1} '.nii.gz' ];
                exec_cmd = ['fslmaths ' temp_seg_L  cur_segname_L] ;
                %Creating a list to be saved:
                obj.Params.tracx_thal2papez.in.seg_list{ii} = obj.Params.tracx_thal2papez.in.temp_list{ii}{1};
                
                
                if exist(cur_segname_L,'file') == 0
                    display([' In proc_tracx2papez(): merging ' cur_segname_L '...'])
                    obj.RunBash(exec_cmd);
                    display('...done');
                end
                %rh:
                cur_segname_R = [ obj.Params.tracx_thal2papez.in.segs_dir 'rh_' ...
                    obj.Params.tracx_thal2papez.in.temp_list{ii}{1} '.nii.gz' ];
                exec_cmd = ['fslmaths ' temp_seg_R cur_segname_R ] ;
                if exist(cur_segname_R) == 0
                    display([' In proc_tracx2papez(): merging ' cur_segname_R '...'])
                    obj.RunBash(exec_cmd);
                    display('...done');
                end
            end
            
            %Creating the txt files:
            %lh:
            obj.Params.tracx_thal2papez.in.lh_txt =  [ obj.Params.tracx_thal2papez.in.segs_dir 'segs_lh.txt' ];
            exec_cmd = [ ' ls -1 ' obj.Params.tracx_thal2papez.in.segs_dir ...
                'lh_* > ' obj.Params.tracx_thal2papez.in.lh_txt  ];
            obj.RunBash(exec_cmd);
            %rh:
            obj.Params.tracx_thal2papez.in.rh_txt =  [ obj.Params.tracx_thal2papez.in.segs_dir 'segs_rh.txt' ];
            exec_cmd = [ ' ls -1 ' obj.Params.tracx_thal2papez.in.segs_dir ...
                'rh_* > ' obj.Params.tracx_thal2papez.in.rh_txt  ];
            obj.RunBash(exec_cmd);
            
            
            
            %Now ready to perform probablistic tractograhy with the
            %parameters of interest:
            %Selecingt our seeds of interes:
            obj.Params.tracx_thal2papez.in.lh_seed = strrep(obj.Params.tracx_thal2papez.in.FSaparc_dir,['_aseg' filesep'],['2009_aseg' filesep 'dwi_fs_Left-Thalamus-Proper.nii.gz']);
            obj.Params.tracx_thal2papez.in.rh_seed = strrep(obj.Params.tracx_thal2papez.in.FSaparc_dir,['_aseg' filesep'],['2009_aseg' filesep 'dwi_fs_Right-Thalamus-Proper.nii.gz']);
            
            %lh:
            exec_cmd = ['probtrackx2 -x ' obj.Params.tracx_thal2papez.in.lh_seed ...
                ' -l --loopcheck -c 0.2 -S 2000 --steplength=0.5 -P 5000' ...
                ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s ' ...
                ' ' obj.Params.tracx_thal2papez.in.bedp_dir filesep 'merged' ...
                ' -m ' obj.Params.tracx_thal2papez.in.bedp_dir filesep 'nodif_brain_mask' ...
                ' --dir=' obj.Params.tracx_thal2papez.in_dir ' --otargetpaths' ...
                ' --targetmasks=' obj.Params.tracx_thal2papez.in.lh_txt ' --os2t'];
            
            
            
            obj.Params.tracx_thal2papez.out.seed2temporal_lh = [ obj.Params.tracx_thal2papez.in_dir 'seeds_to_lh_temporal.nii.gz'] ;
            if exist(obj.Params.tracx_thal2papez.out.seed2temporal_lh) == 0
                obj.RunBash(exec_cmd,44);
            end
            %rh:
            exec_cmd = ['probtrackx2 -x ' obj.Params.tracx_thal2papez.in.rh_seed ...
                ' -l --loopcheck -c 0.2 -S 2000 --steplength=0.5 -P 5000' ...
                ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s ' ...
                ' ' obj.Params.tracx_thal2papez.in.bedp_dir filesep 'merged' ...
                ' -m ' obj.Params.tracx_thal2papez.in.bedp_dir filesep 'nodif_brain_mask' ...
                ' --dir=' obj.Params.tracx_thal2papez.in_dir ' --otargetpaths' ...
                ' --targetmasks=' obj.Params.tracx_thal2papez.in.rh_txt ' --os2t'];
            
            obj.Params.tracx_thal2papez.out.seed2temporal_rh = [ obj.Params.tracx_thal2papez.in_dir 'seeds_to_lh_temporal.nii.gz'] ;
            if exist(obj.Params.tracx_thal2papez.out.seed2temporal_rh,'file') == 0
                obj.RunBash(exec_cmd,44);
            end
            
            %Now find the biggest
            %lh:
            obj.Params.tracx_thal2papez.out.biggest_lh = [ obj.Params.tracx_thal2papez.in_dir 'lh_thal2ctx11.nii.gz'] ;
            if exist(obj.Params.tracx_thal2papez.out.biggest_lh,'file') == 0
                exec_cmd = ['find_the_biggest ' obj.Params.tracx_thal2papez.in_dir  ...
                    'seeds_to_lh_* '  obj.Params.tracx_thal2papez.out.biggest_lh ]
                obj.RunBash(exec_cmd,44);
            end
            %Extracting the values and uploading data
            [check_ok , tmp_lh_vals] = system(['fslstats ' ...
                obj.Params.tracx_thal2papez.out.biggest_lh  ' -h ' num2str([1+numel(obj.Params.tracx_thal2papez.in.seg_list)]) ]);
            if check_ok ~= 0
                error('In  proc_tracx2papez(): failed to parcellate find_the_biggest for thal2ctx11_lh');
            end
            %formating the values for double characters
            floats_tmp_lh_vals=textscan(tmp_lh_vals,'%f');
            obj.Params.tracx_thal2papez.out.biggest_vals_lh = floats_tmp_lh_vals{1};
            
            
            %rh
            obj.Params.tracx_thal2papez.out.biggest_rh = [ obj.Params.tracx_thal2papez.in_dir 'rh_thal2ctx11.nii.gz'] ;
            if exist(obj.Params.tracx_thal2papez.out.biggest_rh,'file') == 0
                exec_cmd = ['find_the_biggest ' obj.Params.tracx_thal2papez.in_dir  ...
                    'seeds_to_rh_* '  obj.Params.tracx_thal2papez.out.biggest_rh ]
                obj.RunBash(exec_cmd,44);
            end
            
            
            %Extracting the values and uploading data
            [check_ok , obj.Params.tracx_thal2papez.out.biggest_vals_rh] = system(['fslstats ' ...
                obj.Params.tracx_thal2papez.out.biggest_rh  ' -h ' num2str([1+numel(obj.Params.tracx_thal2papez.in.seg_list)]) ]);
            if check_ok ~= 0
                error('In  proc_tracx2papez(): failed to parcellate find_the_biggest for thal2ctx11_rh' );
            end
            
            
            %Saving object:
            obj.resave
            
        end
        
        
        
        %%GET DATA Functions
        function obj = getdata_FreeSurfer(obj)
         
            files = {[obj.FS_location filesep obj.sessionname filesep 'stats' filesep 'lh.aparc.stats'];
                [obj.FS_location filesep obj.sessionname filesep 'stats' filesep 'rh.aparc.stats'];
                [obj.FS_location filesep obj.sessionname filesep 'stats' filesep 'aseg.stats']};
            labs = {'lh' 'rh' 'vol'};
            
            
            for zz = 1:3
                tmp_all = ReadInFile(files{zz},'\t',0); tmp = tmp_all(contains('^# ColHeaders',tmp_all):end);
                
                
                
                for ii = 1:50
                    tmp = regexprep(tmp,'  ', ' ');
                end
                
                tmp = regexp(tmp,' ','split');
                T = {};
                for ii = 1:numel(tmp);
                    if isempty(tmp{ii}{end}); tmp{ii} = tmp{ii}(1:end-1); end
                    if ii == 1
                        T(ii,:) = tmp{ii}(3:end);
                    else
                        if isempty(str2num(tmp{ii}{1}))
                            tmp{ii}([2:numel(tmp{ii})]) = num2cell(str2num(char(tmp{ii}([2:numel(tmp{ii})]))));
                        else
                            tmp{ii}([1:4 6:numel(tmp{ii})]) = num2cell(str2num(char(tmp{ii}([1:4 6:numel(tmp{ii})]))));
                        end
                        T(ii,:) = tmp{ii}(1:numel(tmp{ii}));
                    end
                end
                
                %Getting estimated measures 
                if zz == 3 %In aparc.stats
                    col_estimatedMeas=tmp_all(contains('^# Measure',tmp_all));
                    for jj=1:numel(col_estimatedMeas)
                        spl_TICV{jj}=strsplit(col_estimatedMeas{jj},', ');
                        T(end+1,4)=  num2cell(str2num(strrep(spl_TICV{jj}{end-1},',','')));
                        T(end,5) = {(strrep(spl_TICV{jj}{2},'.',''))};
                    end
                end
                
                obj.FSdata.(labs{zz}) = cell2table(T(2:end,:),'VariableName',T(1,:));
            end
            obj.resave;
        end
        
        function obj = getdata_trkland_fx(obj)
            %Get volume data of unclean/cleaned tracts:
            %Volume:
            obj.Trkland.fx.data.lh_unclean_vol = obj.Trkland.Trks.fx_raw_lh.num_uvox;
            obj.Trkland.fx.data.lh_clean_vol = obj.Trkland.Trks.fx_cleantrimmed_lh.num_uvox;
            obj.Trkland.fx.data.lh_trimmedclean_vol = obj.Trkland.Trks.fx_trimmed_lh.num_uvox;
            obj.Trkland.fx.data.lh_cline_HighFA_vol = obj.Trkland.Trks.fx_clinehighFA_lh.num_uvox;
            obj.Trkland.fx.data.lh_cline_HDorff_vol = obj.Trkland.Trks.fx_clineHDorff_lh.num_uvox;
            %METRICS DATA NOW
            %unclean_fx_lh
            obj.Trkland.fx.data.lh_unclean_FA = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,4));
            obj.Trkland.fx.data.lh_unclean_RD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,5));
            obj.Trkland.fx.data.lh_unclean_AxD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,6));
            obj.Trkland.fx.data.lh_unclean_MD = mean(obj.Trkland.Trks.fx_raw_lh.unique_voxels(:,7));
            %Trimmed_clean_hippocing:
            obj.Trkland.fx.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,4));
            obj.Trkland.fx.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,5));
            obj.Trkland.fx.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,6));
            obj.Trkland.fx.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.fx_trimmed_lh.unique_voxels(:,7));
            %Clean_fx:
            obj.Trkland.fx.data.lh_clean_FA = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,4));
            obj.Trkland.fx.data.lh_clean_RD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,5));
            obj.Trkland.fx.data.lh_clean_AxD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,6));
            obj.Trkland.fx.data.lh_clean_MD = mean(obj.Trkland.Trks.fx_cleantrimmed_lh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.fx_clinehighFA_lh,'maxsstrlen')
                obj.Trkland.Trks.fx_clinehighFA_lh = rotrk_centerline(obj.Trkland.Trks.fx_cleantrimmed_lh, 'high_sc','FA');
            end
            obj.Trkland.fx.data.lh_cline_length_highFA=obj.Trkland.Trks.fx_clinehighFA_lh.maxsstrlen;
            obj.Trkland.fx.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,4));
            obj.Trkland.fx.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,5));
            obj.Trkland.fx.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,6));
            obj.Trkland.fx.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_lh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.fx_clineHDorff_lh,'maxsstrlen')
                obj.Trkland.Trks.fx_clineHDorff_lh = rotrk_centerline(obj.Trkland.Trks.fx_cleantrimmed_lh, 'hausdorff');
            end
            obj.Trkland.fx.data.lh_cline_length_HDorff=obj.Trkland.Trks.fx_clineHDorff_lh.maxsstrlen;
            obj.Trkland.fx.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,4));
            obj.Trkland.fx.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,5));
            obj.Trkland.fx.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,6));
            obj.Trkland.fx.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_lh.unique_voxels(:,7));
          
            %RH:
            %Get volume data of unclean/cleaned tracts:
            %Volume:
            obj.Trkland.fx.data.rh_unclean_vol = obj.Trkland.Trks.fx_raw_rh.num_uvox;
            obj.Trkland.fx.data.rh_clean_vol = obj.Trkland.Trks.fx_cleantrimmed_rh.num_uvox;
            obj.Trkland.fx.data.rh_trimmedclean_vol = obj.Trkland.Trks.fx_trimmed_rh.num_uvox;
            obj.Trkland.fx.data.rh_cline_HighFA_vol = obj.Trkland.Trks.fx_clinehighFA_rh.num_uvox;
            obj.Trkland.fx.data.rh_cline_HDorff_vol = obj.Trkland.Trks.fx_clineHDorff_rh.num_uvox;
            
            %METRICS DATA NOW
            %unclean_fx_rh
            obj.Trkland.fx.data.rh_unclean_FA = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,4));
            obj.Trkland.fx.data.rh_unclean_RD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,5));
            obj.Trkland.fx.data.rh_unclean_AxD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,6));
            obj.Trkland.fx.data.rh_unclean_MD = mean(obj.Trkland.Trks.fx_raw_rh.unique_voxels(:,7));
            %Trimmed_clean_hippocing:
            obj.Trkland.fx.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,4));
            obj.Trkland.fx.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,5));
            obj.Trkland.fx.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,6));
            obj.Trkland.fx.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.fx_trimmed_rh.unique_voxels(:,7));
            %Clean_fx:
            obj.Trkland.fx.data.rh_clean_FA = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,4));
            obj.Trkland.fx.data.rh_clean_RD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,5));
            obj.Trkland.fx.data.rh_clean_AxD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,6));
            obj.Trkland.fx.data.rh_clean_MD = mean(obj.Trkland.Trks.fx_cleantrimmed_rh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.fx_clinehighFA_rh,'maxsstrlen')
                obj.Trkland.Trks.fx_clinehighFA_rh = rotrk_centerline(obj.Trkland.Trks.fx_cleantrimmed_rh, 'high_sc','FA');
            end
            obj.Trkland.fx.data.rh_cline_length_HDorff=obj.Trkland.Trks.fx_clinehighFA_rh.maxsstrlen;
            obj.Trkland.fx.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,4));
            obj.Trkland.fx.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,5));
            obj.Trkland.fx.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,6));
            obj.Trkland.fx.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.fx_clinehighFA_rh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.fx_clineHDorff_rh,'maxsstrlen')
                obj.Trkland.Trks.fx_clineHDorff_rh = rotrk_centerline(obj.Trkland.Trks.fx_cleantrimmed_rh, 'hausdorff');
            end
            obj.Trkland.fx.data.rh_cline_length_HDorff=obj.Trkland.Trks.fx_clineHDorff_rh.maxsstrlen;
            obj.Trkland.fx.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,4));
            obj.Trkland.fx.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,5));
            obj.Trkland.fx.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,6));
            obj.Trkland.fx.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.fx_clineHDorff_rh.unique_voxels(:,7));
            
            
            
            %Saving to *.mat:
            obj.resave();
        end
        
        function obj = getdata_trkland_hippocing(obj)
            %METRICS DATA NOW
            %LH:
            %unclean_hippocing_lh
            obj.Trkland.hippocing.data.lh_unclean_FA = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.lh_unclean_RD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.lh_unclean_AxD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.lh_unclean_MD = mean(obj.Trkland.Trks.raw_hippocing_lh.unique_voxels(:,7));
            %Clean_hippocing:
            obj.Trkland.hippocing.data.lh_clean_FA = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.lh_clean_RD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.lh_clean_AxD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.lh_clean_MD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_lh.unique_voxels(:,7));
            %Trimmed_clean_hippocing:
            obj.Trkland.hippocing.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.hippocing_trimmed_lh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.hippocing_clinehighFA_lh,'maxsstrlen')
                obj.Trkland.Trks.hippocing_clinehighFA_lh= rotrk_centerline(obj.Trkland.Trks.hippocing_cleantrimmed_lh, 'high_sc','FA');
            end
            obj.Trkland.hippocing.data.lh_cline_length_highFA=obj.Trkland.Trks.hippocing_clinehighFA_lh.maxsstrlen;
            obj.Trkland.hippocing.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_lh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.hippocing_clineHDorff_lh,'maxsstrlen')
                obj.Trkland.Trks.hippocing_clineHDorff_lh = rotrk_centerline(obj.Trkland.Trks.hippocing_cleantrimmed_lh, 'hausdorff');
            end
            obj.Trkland.hippocing.data.lh_cline_length_HDorff=obj.Trkland.Trks.hippocing_clineHDorff_lh.maxsstrlen;
            obj.Trkland.hippocing.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_lh.unique_voxels(:,7));
            
            
            %RH:
            %Get volume data of unclean/cleaned tracts:
            %Volume:
            obj.Trkland.hippocing.data.rh_unclean_vol = obj.Trkland.Trks.raw_hippocing_rh.num_uvox;
            obj.Trkland.hippocing.data.rh_clean_vol = obj.Trkland.Trks.hippocing_cleantrimmed_rh.num_uvox;
            obj.Trkland.hippocing.data.rh_trimmedclean_vol = obj.Trkland.Trks.hippocing_trimmed_rh.num_uvox;
            obj.Trkland.hippocing.data.rh_cline_HighFA_vol = obj.Trkland.Trks.hippocing_clinehighFA_rh.num_uvox;
            obj.Trkland.hippocing.data.rh_cline_HDorff_vol = obj.Trkland.Trks.hippocing_clineHDorff_rh.num_uvox;
            %METRICS DATA NOW
            %unclean_hippocing_rh
            obj.Trkland.hippocing.data.rh_unclean_FA = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.rh_unclean_RD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.rh_unclean_AxD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.rh_unclean_MD = mean(obj.Trkland.Trks.raw_hippocing_rh.unique_voxels(:,7));
            %Clean_hippocing:
            obj.Trkland.hippocing.data.rh_clean_FA = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.rh_clean_RD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.rh_clean_AxD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.rh_clean_MD = mean(obj.Trkland.Trks.hippocing_cleantrimmed_rh.unique_voxels(:,7));
            %Trimmed_clean_hippocing:
            obj.Trkland.hippocing.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.hippocing_trimmed_rh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.hippocing_clinehighFA_rh,'maxsstrlen')
                obj.Trkland.Trks.hippocing_clinehighFA_rh= rotrk_centerline(obj.Trkland.Trks.hippocing_cleantrimmed_rh, 'high_sc','FA');
            end
            obj.Trkland.hippocing.data.rh_cline_length_highFA=obj.Trkland.Trks.hippocing_clinehighFA_rh.maxsstrlen;
            obj.Trkland.hippocing.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.hippocing_clinehighFA_rh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.hippocing_clineHDorff_rh,'maxsstrlen')
                obj.Trkland.Trks.hippocing_clineHDorff_rh = rotrk_centerline(obj.Trkland.Trks.hippocing_cleantrimmed_rh, 'hausdorff');
            end
            obj.Trkland.hippocing.data.rh_cline_length_HDorff=obj.Trkland.Trks.hippocing_clineHDorff_rh.maxsstrlen;
            obj.Trkland.hippocing.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,4));
            obj.Trkland.hippocing.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,5));
            obj.Trkland.hippocing.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,6));
            obj.Trkland.hippocing.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.hippocing_clineHDorff_rh.unique_voxels(:,7));
            
            
            %Resave:
            obj.resave();
            
        end
        
        function obj = getdata_trkland_cingulum(obj)
            %Get volume data of unclean/cleaned tracts:
            %Volume:
            obj.Trkland.cingulum.data.lh_unclean_vol = obj.Trkland.Trks.raw_cingulum_lh.num_uvox;
            obj.Trkland.cingulum.data.lh_trimmedclean_vol = obj.Trkland.Trks.cingulum_trimmed_lh.num_uvox;
            obj.Trkland.cingulum.data.lh_clean_vol = obj.Trkland.Trks.cingulum_cleantrimmed_lh.num_uvox;
            obj.Trkland.cingulum.data.lh_cline_HighFA_vol = obj.Trkland.Trks.cingulum_clinehighFA_lh.num_uvox;
            obj.Trkland.cingulum.data.lh_cline_HDorff_vol = obj.Trkland.Trks.cingulum_clineHDorff_lh.num_uvox;
            %METRICS DATA NOW
            %unclean_cingulum_lh
            obj.Trkland.cingulum.data.lh_unclean_FA = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.lh_unclean_RD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.lh_unclean_AxD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.lh_unclean_MD = mean(obj.Trkland.Trks.raw_cingulum_lh.unique_voxels(:,7));
            %Clean_cingulum:
            obj.Trkland.cingulum.data.lh_clean_FA = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.lh_clean_RD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.lh_clean_AxD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.lh_clean_MD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_lh.unique_voxels(:,7));
            %Trimmed_clean_cingulum:
            obj.Trkland.cingulum.data.lh_trimmedclean_FA = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.lh_trimmedclean_RD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.lh_trimmedclean_AxD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.lh_trimmedclean_MD = mean(obj.Trkland.Trks.cingulum_trimmed_lh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.cingulum_clinehighFA_lh,'maxsstrlen')
                obj.Trkland.Trks.cingulum_clinehighFA_lh= rotrk_centerline(obj.Trkland.Trks.cingulum_cleantrimmed_lh, 'high_sc','FA');
            end
            obj.Trkland.cingulum.data.lh_cline_length_highFA=obj.Trkland.Trks.cingulum_clinehighFA_lh.maxsstrlen;
            obj.Trkland.cingulum.data.lh_cline_FA_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.lh_cline_RD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.lh_cline_AxD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.lh_cline_MD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_lh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.cingulum_clineHDorff_lh,'maxsstrlen')
                obj.Trkland.Trks.cingulum_clineHDorff_lh = rotrk_centerline(obj.Trkland.Trks.cingulum_cleantrimmed_lh, 'hausdorff');
            end
            obj.Trkland.cingulum.data.lh_cline_length_HDorff=obj.Trkland.Trks.cingulum_clineHDorff_lh.maxsstrlen;
            obj.Trkland.cingulum.data.lh_cline_FA_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.lh_cline_RD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.lh_cline_AxD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.lh_cline_MD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_lh.unique_voxels(:,7));
            
            
            %RH:
            %Get volume data of unclean/cleaned tracts:
            %Volume:
            obj.Trkland.cingulum.data.rh_unclean_vol = obj.Trkland.Trks.raw_cingulum_rh.num_uvox;
            obj.Trkland.cingulum.data.rh_clean_vol = obj.Trkland.Trks.cingulum_cleantrimmed_rh.num_uvox;
            obj.Trkland.cingulum.data.rh_trimmedclean_vol = obj.Trkland.Trks.cingulum_trimmed_rh.num_uvox;
            obj.Trkland.cingulum.data.rh_cline_HighFA_vol = obj.Trkland.Trks.cingulum_clinehighFA_rh.num_uvox;
            obj.Trkland.cingulum.data.rh_cline_HDorff_vol = obj.Trkland.Trks.cingulum_clineHDorff_rh.num_uvox;
            %METRICS DATA NOW
            %unclean_cingulum_rh
            obj.Trkland.cingulum.data.rh_unclean_FA = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.rh_unclean_RD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.rh_unclean_AxD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.rh_unclean_MD = mean(obj.Trkland.Trks.raw_cingulum_rh.unique_voxels(:,7));
            %Clean_cingulum:
            obj.Trkland.cingulum.data.rh_clean_FA = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.rh_clean_RD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.rh_clean_AxD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.rh_clean_MD = mean(obj.Trkland.Trks.cingulum_cleantrimmed_rh.unique_voxels(:,7));
            %Trimmed_clean_cingulum:
            obj.Trkland.cingulum.data.rh_trimmedclean_FA = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.rh_trimmedclean_RD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.rh_trimmedclean_AxD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.rh_trimmedclean_MD = mean(obj.Trkland.Trks.cingulum_trimmed_rh.unique_voxels(:,7));
            %Cline_HighFA
            if ~isfield(obj.Trkland.Trks.cingulum_clinehighFA_rh,'maxsstrlen')
                obj.Trkland.Trks.cingulum_clinehighFA_rh= rotrk_centerline(obj.Trkland.Trks.cingulum_cleantrimmed_rh, 'high_sc','FA');
            end
            obj.Trkland.cingulum.data.rh_cline_length_highFA=obj.Trkland.Trks.cingulum_clinehighFA_rh.maxsstrlen;
            obj.Trkland.cingulum.data.rh_cline_FA_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.rh_cline_RD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.rh_cline_AxD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.rh_cline_MD_highFA = mean(obj.Trkland.Trks.cingulum_clinehighFA_rh.unique_voxels(:,7));
            %Cline_HDorff
            if ~isfield(obj.Trkland.Trks.cingulum_clineHDorff_rh,'maxsstrlen')
                obj.Trkland.Trks.cingulum_clineHDorff_rh = rotrk_centerline(obj.Trkland.Trks.cingulum_cleantrimmed_rh, 'hausdorff');
            end
            obj.Trkland.cingulum.data.rh_cline_length_HDorff=obj.Trkland.Trks.cingulum_clineHDorff_rh.maxsstrlen;
            obj.Trkland.cingulum.data.rh_cline_FA_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,4));
            obj.Trkland.cingulum.data.rh_cline_RD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,5));
            obj.Trkland.cingulum.data.rh_cline_AxD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,6));
            obj.Trkland.cingulum.data.rh_cline_MD_HDorff = mean(obj.Trkland.Trks.cingulum_clineHDorff_rh.unique_voxels(:,7));
            
            %RESAVING:
            obj.resave();
            
            
        end
        
        %%%%%%%%%%%%%%%%%%% END Data Post-Processing Methods %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% END of Post-Processing Methods %%%%%%%%%%%%%%%%%
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
        
        function obj = RunBash(obj,exec_cmd, exit_status)
            %Code values:
            %   44  --> Show output
            %   144 --> Show output and exist_status is 1
            
            
            if nargin < 3 %This will allow us to pass other exit status, such as 1 for DSISTUDIO in GQI
                exit_status = 0 ;
            end
            
            if exit_status == 44 || exit_status == 144 % 44 codes for non-comment when running system:
                [ sys_success , sys_error ] = system(exec_cmd,'-echo') ;
                if exit_status == 44
                    exit_status = 0; %return exit_status to zero.
                else
                    exit_status = 1; %return exit_status to zero.
                end
            else
                [ sys_success , sys_error ] = system(exec_cmd) ;
            end
            
            if sys_success==exit_status; disp('');
            else
                fprintf('\n\n===================================');
                fprintf('===================================\n\n');
                fprintf(['\n \t\tError in exec_cmd: \n ' exec_cmd '\n\n\n' ]);
                fprintf(['\n \t\tError output by system(exec_cmd) is:\n ' sys_error '\n'])
                fprintf('\n\n===================================');
                fprintf('===================================\n\n');
                error('Stopping now....');
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
                    system(['mkdir -p ' obj.root ]);
                catch
                    disp([ 'Trying to create /DWIs/ in:' obj.root ...
                        ' Maybe some permission issues?'])
                end
            end
        end
        
        function obj = RefreshFields(obj,whatParam,direction)
            if nargin <3
                direction = 'bil' ; %assuming bilateral directionality if not 3rd argument given
            end
            switch whatParam
                case 'hippocing'
                    fields_out = fields(obj.Trkland.hippocing.out);
                    fields_data = fields(obj.Trkland.hippocing.data);
                    for ii=1:numel(fields_out)
                        if strcmp(direction,'bil')
                            obj.Trkland.hippocing.out.(fields_out{ii}) = '' ;
                        else
                            if ~isempty(strfind(fields_out{ii},direction))
                                obj.Trkland.hippocing.out.(fields_out{ii}) = '' ;
                            end
                        end
                    end
                    for ii=1:numel(fields_data)
                        if strcmp(direction,'bil')
                            obj.Trkland.hippocing.data.(fields_data{ii}) = [];
                        else
                            if ~isempty(strfind(fields_data{ii},direction))
                                obj.Trkland.hippocing.data(fields_data{ii}) = [] ;
                            end
                        end
                    end
                case 'atr'
                    fields_out = fields(obj.Trkland.atr.out);
                    fields_data = fields(obj.Trkland.atr.data);
                    for ii=1:numel(fields_out)
                        if strcmp(direction,'bil')
                            obj.Trkland.atr.out.(fields_out{ii}) = '' ;
                        else
                            if ~isempty(strfind(fields_out{ii},direction))
                                obj.Trkland.atr.out.(fields_out{ii}) = '' ;
                            end
                        end
                    end
                    for ii=1:numel(fields_data)
                        if strcmp(direction,'bil')
                            obj.Trkland.atr.data.(fields_data{ii}) = [];
                        else
                            if ~isempty(strfind(fields_data{ii},direction))
                                obj.Trkland.atr.data(fields_data{ii}) = [] ;
                            end
                        end
                    end
                case 'fx'
                    fields_out = fields(obj.Trkland.fx.out);
                    fields_data = fields(obj.Trkland.fx.data);
                    for ii=1:numel(fields_out)
                        if strcmp(direction,'bil')
                            obj.Trkland.fx.out.(fields_out{ii}) = '' ;
                        else
                            if ~isempty(strfind(fields_out{ii},direction))
                                obj.Trkland.fx.out.(fields_out{ii}) = '' ;
                            end
                        end
                    end
                    for ii=1:numel(fields_data)
                        if strcmp(direction,'bil')
                            obj.Trkland.fx.data.(fields_data{ii}) = [];
                        else
                            if ~isempty(strfind(fields_data{ii},direction))
                                obj.Trkland.fx.data.(fields_data{ii}) = [] ;
                            end
                        end
                    end
                case 'cingulum'
                    fields_out = fields(obj.Trkland.cingulum.out);
                    fields_data = fields(obj.Trkland.cingulum.data);
                    for ii=1:numel(fields_out)
                        if strcmp(direction,'bil')
                            obj.Trkland.cingulum.out.(fields_out{ii}) = '' ;
                        else
                            if ~isempty(strfind(fields_out{ii},direction))
                                obj.Trkland.cingulum.out.(fields_out{ii}) = '' ;
                            end
                        end
                    end
                    for ii=1:numel(fields_data)
                        if strcmp(direction,'bil')
                            obj.Trkland.cingulum.data.(fields_data{ii}) = [];
                        else
                            if ~isempty(strfind(fields_data{ii},direction))
                                obj.Trkland.cingulum.data.(fields_data{ii}) = [] ;
                            end
                        end
                    end
                otherwise
                    error(['In RefreshFields(). I do not understand the argument: --> ' whatParam ' <-- please check']);
            end
            obj.resave;
        end
        function obj = UploadData_DWI(obj)
            id = obj.sessionname;
            if isempty(id);
                disp('No Session_ID.  Cannot upload data');
                return
            end
            
            if isnumeric(id)
                id = num2str(id);
            end
            
            %%Select current SessionID
            dctl_cmd = [ 'SELECT MRI_Session_ID FROM Sessions.MRI  WHERE ' ' MRI_Session_Name = ''' id '''' ];
            cur_DC_ID = DataCentral(dctl_cmd);
            
            %%Eddymotion uploading
            MOTION_fields=fields(obj.Params.EddyMotion.out.vals);
            for ii=1:numel(MOTION_fields)
                if size(obj.Params.EddyMotion.out.vals.(MOTION_fields{ii}){1},2) == 9 %9 is the number of characters for this double type variable
                    disp(''); %Skip uploading because it already exists!
                else
                    fprintf(['\nUploading motion value: ' MOTION_fields{ii} '(' ...
                        strtrim(cell2char(obj.Params.EddyMotion.out.vals.(MOTION_fields{ii})))  ') for ' obj.sessionname ]);
                    dctl_cmd = [ ' SELECT MRI_skelDWI_motion_eddyres_' MOTION_fields{ii} ...
                        ' FROM MRI.skelDWI  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                    check_dctl_cmd = DataCentral(dctl_cmd);
                    if isempty(check_dctl_cmd.(['MRI_skelDWI_motion_eddyres_' MOTION_fields{ii}]))
                        fprintf(['Motion values: ' MOTION_fields{ii} ' is: '   ])
                        dctl_cmd = [ 'INSERT INTO MRI.skelDWI (MRI_Session_ID,  MRI_skelDWI_motion_eddyres_' MOTION_fields{ii} ') ' ...
                            ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' strtrim(cell2char(obj.Params.EddyMotion.out.vals.(MOTION_fields{ii}))) ')'   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    elseif isnan(check_dctl_cmd.(['MRI_skelDWI_motion_eddyres_' MOTION_fields{ii}]))
                        fprintf(['Skel TOI ==> Uploading to DataCentral: ' id ' and TOI: ' MOTION_fields{ii}  ])
                        dctl_cmd = [ 'UPDATE MRI.skelDWI SET MRI_skelDWI_motion_eddyres_' MOTION_fields{ii} ...
                            ' = ''' strtrim(cell2char(obj.Params.EddyMotion.out.vals.(MOTION_fields{ii}))) ''' WHERE MRI_Session_ID =  ' ...
                            num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    end
                end
            end
            fprintf('...done\n');
            
            %%Skel Values uploading
            TOI_fields=fields(obj.Params.Skel_TOI.out);
            fprintf('\n Uploading: ');
            
            for ii=1:numel(TOI_fields)
                dctl_cmd = [ ' SELECT MRI_skelDWI_' TOI_fields{ii} ...
                    ' FROM MRI.skelDWI  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                %  dctl_cmd = [ ' SELECT MRI_skelDWI_' TOI_fields{ii} ...
                %       ' FROM rdp20.DWI_TBSS_SKEL_VALS  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                
                check_dctl_cmd = DataCentral(dctl_cmd);
                
                %Add 1000 to AxD, RD, and MD
                if isempty(strfind(TOI_fields{ii},'FA'))
                    cur_value=strtrim(num2str(str2num(obj.Params.Skel_TOI.out.(TOI_fields{ii}))*1000));
                else
                    cur_value=strtrim(obj.Params.Skel_TOI.out.(TOI_fields{ii}));
                end
                
                if isempty(check_dctl_cmd.(['MRI_skelDWI_' TOI_fields{ii}]))
                    fprintf(['Skel TOI ==> Uploading to DataCentral: ' id ' and TOI: ' TOI_fields{ii}  ])
                    %       dctl_cmd = [ 'INSERT INTO rdp20.DWI_TBSS_SKEL_VALS (MRI_Session_ID,  MRI_skelDWI_' TOI_fields{ii} ') ' ...
                    %          ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' strtrim(obj.Params.Skel_TOI.out.(TOI_fields{ii})) ')'   ] ;
                    dctl_cmd = [ 'INSERT INTO MRI.skelDWI (MRI_Session_ID,  MRI_skelDWI_' TOI_fields{ii} ') ' ...
                        ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' cur_value ')'   ] ;
                    DataCentral(dctl_cmd);
                    fprintf('...done\n');
                    %elseif isnan(check_dctl_cmd.(['MRI_skelDWI_' TOI_fields{ii}]))
                else
                    %fprintf(['Skel TOI ==> Uploading to DataCentral: ' id ' and TOI: ' TOI_fields{ii}  ])
                    fprintf([TOI_fields{ii} '...' ])
                    if mod(ii,4) == 0
                        fprintf('\n');
                    end
                    %    dctl_cmd = [ 'UPDATE rdp20.DWI_TBSS_SKEL_VALS SET MRI_skelDWI_' TOI_fields{ii} ...
                    %    ' = ''' strtrim(obj.Params.Skel_TOI.out.(TOI_fields{ii})) ''' WHERE MRI_Session_ID =  ' ...
                    %    num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                    dctl_cmd = [ 'UPDATE MRI.skelDWI SET MRI_skelDWI_' TOI_fields{ii} ...
                        ' = ''' cur_value ''' WHERE MRI_Session_ID =  ' ...
                        num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                    DataCentral(dctl_cmd);
                    %fprintf('...done\n');
                end
                % DataCentral(['UPDATE PIB.noT1_SUVR_new SET ' Q ' WHERE PIB_Session_ID = "' id '"']);
                
                %%% upload ROI data
                obj.resave;
            end
            disp('skelTOIs values have been uploaded to DataCentral');
        end
        
        
        
        function obj = UploadData_Trkland(obj)
            id = obj.sessionname;
            if isempty(id);
                disp('No Session_ID.  Cannot upload data');
                return
            end
            
            if isnumeric(id)
                id = num2str(id);
            end
            
            %%Select current SessionID
            dctl_cmd = [ 'SELECT MRI_Session_ID FROM Sessions.MRI  WHERE ' ' MRI_Session_Name = ''' id '''' ];
            cur_DC_ID = DataCentral(dctl_cmd);
            
            %%Creating all the files that will be uploaded:
            FX_fields = fields(obj.Trkland.fx.data);
            CING_fields = fields(obj.Trkland.cingulum.data);
            HIPPOCING_fields=fields(obj.Trkland.hippocing.data);
            
            %STARTING THE UPLOADING SEQUENCE:
            
            fprintf('\n Uploading: ');
            
            %FOR FX:
            for ii=1:numel(FX_fields)
                dctl_cmd = [ ' SELECT fx_' FX_fields{ii} ...
                    ' FROM rdp20.TRKLAND  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                check_dctl_cmd = DataCentral(dctl_cmd);
                
                %Add 1000 to AxD, RD, and MD
                if ~isempty(obj.Trkland.fx.data.(FX_fields{ii}))
                    if ~isempty(strfind(FX_fields{ii},'RD')) || ~isempty(strfind(FX_fields{ii},'MD')) || ~isempty(strfind(FX_fields{ii},'AxD'))
                        cur_value=num2str(obj.Trkland.fx.data.(FX_fields{ii})*1000);
                    else
                        cur_value=num2str(obj.Trkland.fx.data.(FX_fields{ii}));
                    end
                    
                    if isempty(check_dctl_cmd.(['fx_' FX_fields{ii}]))
                        fprintf(['FX Trkland ==> Inserting to DataCentral: ' id ' on: ' FX_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'INSERT INTO rdp20.TRKLAND (MRI_Session_ID,  fx_' FX_fields{ii} ') ' ...
                            ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' cur_value ')'   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    else
                        fprintf(['FX Trkland ==> Updating to DataCentral: ' id ' on: ' FX_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'UPDATE rdp20.TRKLAND SET fx_' FX_fields{ii} ...
                            ' = ''' cur_value ''' WHERE MRI_Session_ID =  ' ...
                            num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    end
                end
                %%% upload ROI data
                obj.resave;
            end
            
            
            %FOR CING:
            for ii=1:numel(CING_fields)
                dctl_cmd = [ ' SELECT cing_' CING_fields{ii} ...
                    ' FROM rdp20.TRKLAND  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                check_dctl_cmd = DataCentral(dctl_cmd);
                display(num2str(ii))
                if ~isempty(obj.Trkland.cingulum.data.(CING_fields{ii}))
                    %Add 1000 to AxD, RD, and MD
                    if ~isempty(strfind(CING_fields{ii},'RD')) || ~isempty(strfind(CING_fields{ii},'MD')) || ~isempty(strfind(CING_fields{ii},'AxD'))
                        cur_value=num2str(obj.Trkland.cingulum.data.(CING_fields{ii})*1000);
                    else
                        cur_value=num2str(obj.Trkland.cingulum.data.(CING_fields{ii}));
                    end
                    
                    if isempty(check_dctl_cmd.(['cing_' CING_fields{ii}]))
                        fprintf(['CINGULUM Trkland ==> Inserting to DataCentral: ' id ' on: ' CING_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'INSERT INTO rdp20.TRKLAND (MRI_Session_ID,  cing_' CING_fields{ii} ') ' ...
                            ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' cur_value ')'   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    else
                        fprintf(['CINGULUM Trkland ==> Updating to DataCentral: ' id ' on: ' CING_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'UPDATE rdp20.TRKLAND SET cing_' CING_fields{ii} ...
                            ' = ''' cur_value ''' WHERE MRI_Session_ID =  ' ...
                            num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    end
                end
                %%% upload ROI data
                obj.resave;
            end
            
            %FOR HIPPOCING:
            for ii=1:numel(HIPPOCING_fields)
                dctl_cmd = [ ' SELECT hippocing_' HIPPOCING_fields{ii} ...
                    ' FROM rdp20.TRKLAND  WHERE MRI_Session_ID = ' num2str(cur_DC_ID.MRI_Session_ID)  ];
                check_dctl_cmd = DataCentral(dctl_cmd);
                display(num2str(ii))
                if ~isempty(obj.Trkland.hippocing.data.(HIPPOCING_fields{ii}))
                    %Add 1000 to AxD, RD, and MD
                    if ~isempty(strfind(HIPPOCING_fields{ii},'RD')) || ~isempty(strfind(HIPPOCING_fields{ii},'MD')) || ~isempty(strfind(HIPPOCING_fields{ii},'AxD'))
                        cur_value=num2str(obj.Trkland.hippocing.data.(HIPPOCING_fields{ii})*1000);
                    else
                        cur_value=num2str(obj.Trkland.hippocing.data.(HIPPOCING_fields{ii}));
                    end
                    
                    if isempty(check_dctl_cmd.(['hippocing_' HIPPOCING_fields{ii}]))
                        fprintf(['HIPPOCINGULUM Trkland ==> Inserting to DataCentral: ' id ' on: ' HIPPOCING_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'INSERT INTO rdp20.TRKLAND (MRI_Session_ID,  hippocing_' HIPPOCING_fields{ii} ') ' ...
                            ' values ( ' num2str(cur_DC_ID.MRI_Session_ID) ',' cur_value ')'   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    else
                        fprintf(['HIPPOCINGULUM Trkland ==> Updating to DataCentral: ' id ' on: ' HIPPOCING_fields{ii} '=' cur_value  ])
                        dctl_cmd = [ 'UPDATE rdp20.TRKLAND SET hippocing_' HIPPOCING_fields{ii} ...
                            ' = ''' cur_value ''' WHERE MRI_Session_ID =  ' ...
                            num2str(cur_DC_ID.MRI_Session_ID)   ] ;
                        DataCentral(dctl_cmd);
                        fprintf('...done\n');
                    end
                end
                %%% upload ROI data
                obj.resave;
            end
            disp('skelTOIs values have been uploaded to DataCentral');
        end
        
        %Methods that can be accessed by other methods only
        %%NOT TESTED --> proc_as_coreg()
        function obj = proc_as_coreg(obj,source,target,list_toxform,coreg_tech)
            %              fprintf('\n%s\n', 'PERFORMING COREGISTRATION');
            %              wasRun = false;
            %
            %              if nargin<2
            %                  error(['Please specify a source (movable image) when ' ...
            %                      'running proc_as_coreg. One arguments missing']);
            %              end
            %
            %              if iscell(source)
            %                  source = char(source);
            %              end
            %
            %              if nargin<3
            %                  error(['Please specify a target (reference image) when ' ...
            %                      'running proc_as_coreg. Two argument missing']);
            %              end
            %
            %             if nargin<4
            %                  fprintf('Coregistering only source image, not list_toxform provided...')
            %                  list_toxform = '';
            %              end
            %
            %              if nargin<5
            %                  coreg_tech = 'spm';
            %              end
            %
            %              if iscell(target)
            %                  target = char(target);
            %              end
            %
            %              [a, ~, ~] = fileparts(source);
            %              outpath1 = obj.getPath(a,source);
            %
            %
            %              switch coreg_tech
            %                  case 'spm'
            %                      disp('USING SPM:');
            %                      if exist([outpath2 'CoReg.mat'],'file')==0
            %                          wasRun=true;
            %
            %                          VG = spm_vol(target);
            %                          VF = spm_vol(source);
            %
            %                          x = spm_coreg(VG,VF,obj.Params.Coreg.in.spm);
            %                          M  = spm_matrix(x);
            %                          save([outpath2 'CoReg.mat'],'x','M');
            %                      else
            %                          load([outpath2 'CoReg.mat']);
            %                      end
            %                      obj.Params.Coreg.out.regfile = [outpath2 'CoReg.mat'];
            %                      disp('Coregistration Estimation has been performed:');
            %
            %                      [a1 b1 c1] = fileparts(source);
            %                      new = [outpath1  obj.Params.Coreg.in.prefix b1 c1];
            %                      if exist(new,'file')==0
            %                          wasRun=true;
            %                          copyfile(source,new);
            %                          MM = spm_get_space(new);
            %                          spm_get_space(new, M\MM);
            %                          if obj.Params.Coreg.in.reslice
            %                              h1 = spm_vol(obj.Params.Coreg.in.target);
            %                              h2 = spm_vol(new);
            %                              m = resizeVol2(h2,h1,obj.Params.Coreg.in.resampopt);
            %                              [a b c] = fileparts(new);
            %                              h1.fname = [a filesep 'rs_' b c];
            %                              h1.dt = h2.dt;
            %                              spm_write_vol(h1,m);
            %                          end
            %                      end
            %                      obj.Params.Coreg.out.regimage = new;
            %                      disp('Coregistration has been applied to the mean image:');
            %
            %                  case 'bbreg'
            %                      error('Not implemented yet...');
            % %                      disp('    USING BBREGISTER:');
            % %                      [a1 b1 c1] = fileparts(source);
            % %
            % %                      reg = [outpath2 'reg_' b1 '.dat'];
            % %                      if exist(reg,'file')==0
            % %                          wasRun = true;
            % %                          cmd = ['bbregister --s  '  obj.fsdir ' --mov ' source ' --init-' obj.Params.FS_Params.bbreg ' --' obj.Params.FS_Params.conweight ' --reg ' reg];
            % %                          runFS(cmd,obj.fsdir);
            % %                      end
            % %                      disp('Coregistration Estimation has been performed:');
            % %                      obj.Params.Coreg.out.regfile = reg;
            % %
            % %                      new = [outpath1 obj.Params.Coreg.in.prefix b1 c1];
            % %                      obj.Params.Coreg.out.regimage = new;
            % %
            % %                      if exist(new,'file')==0
            % %                          wasRun = true;
            % %                          if obj.Params.Coreg.in.reslice
            % %                              cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' source ' --reg ' reg ' --fstarg --o ' new];
            % %                          else
            % %                              cmd = ['mri_vol2vol --s ' obj.fsubj ' --mov ' source ' --reg ' reg ' --fstarg --no-resample --o ' new];
            % %                          end
            % %                          runFS(cmd,obj.fsdir);
            % %                      end
            % %
            % %                      disp('coregistration with bbregister is complete');
            %                  case 'iuw'
            %                      disp('Perform this step with the IUW module');
            %                      return;
            %                  otherwise
            %                      disp('unknown coregistration style option');
            %                      return
            %            end
        end
        
        
        function obj = remove_trkland_fields(obj,curTRK)
            temp_fields = fieldnames(obj.Trkland.Trks.(curTRK));
            for ii=1:numel(temp_fields)
                display(temp_fields{ii});
                if ( strcmp(temp_fields{ii},'filename') || strcmp(temp_fields{ii},'trk_name') ) || strcmp(temp_fields{ii},'id')
                    donothing=1;
                else
                    obj.Trkland.Trks.(curTRK).(temp_fields{ii}) = [];
                end
            end
            clear donothing;
        end
        
        
        function obj = proc_apply_reservsenorm_new(obj,fn,regfile,targ)
            %        obj.Params.ApplyReverseNormNew.in.movefiles = '/autofs/eris/bang/ADRC/Sessions/141015_8CS00118/restingState/';
            % obj.Params.ApplyReverseNormNew.in.fn = '/cluster/brutha/MATLAB_Scripts/PET/FROIs/DMN_PCC.nii';
            % obj.Params.ApplyReverseNormNew.in.targ = '/autofs/eris/bang/ADRC/FreeSurfer6.0/141015_8CS00118/mri/spm12/orig.nii';
            % obj.Params.ApplyReverseNormNew.in.regfile = obj.Params.spmT1_Proc.out.iregfile;
            %
            %
            % obj.proc_apply_reservsenorm_new;
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
            
            if isempty(obj.Params.ApplyReverseNormNew.in.regfile)
                regfile = ''
            end
            
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
        
        
    end
end


