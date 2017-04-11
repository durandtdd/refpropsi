function varargout = refpropsi(prop_req, prop1, value1, prop2, value2, fluid, x, out_type, warning_flag)
%REFPROPSI Compute thermodynamic properties with Refprop
%   Parameters
%   ----------
%       prop_req: str
%           List of properties required
%       prop1: str
%           Property 1 name
%       value1: double
%           Value for property 1
%       prop2: str
%           Property 2 name
%       value2: double
%           Value for property 2
%       fluid: str or cell of str
%           Fluid composition name(s)
%       x: array [optional, default=1]
%           Fluid compostion
%       out_type: str [optional, default=''] {'array','cell','struct',''}
%           Output type (see examples)
%       warning_flag: bool [optional, default=False]
%           If true, show refprop warnings
%     
%
%   Returns
%   -------
%       varargout: Variables or array or cell or struct
%           Requested properties values
%
%
%   Properties
%   ----------
%       In
%       --
%       C   Critical point                  [null]
%
%       In/out:
%       -------
%       D   Density                         [kg/m**3]
%       U   Internal energy                 [J/kg]
%       H   Enthalpy                        [J/kg]
%       P   Pressure                        [Pa]
%       Q   Gas mass fraction               [null]
%       S   Entropy                         [J/kg/K]
%       T   Temperature                     [K]
%
%       Out
%       ---
%       A   Speed of sound                  [m/s]
%       C   Heat capacity cp                [J/kg/K]
%       I   Surface tension                 [N/m]
%       K   Heat capacities ratio (cp/cv)   [J/kg/K]
%       L   Thermal conductivity            [W/m/K]
%       M   Nolar mass                      [kg/mol]
%       O   Heat capacity cv                [J/kg/K]
%       V   Viscosity                       [Pa.s]
%       Z   Compressibilty factor           [null]
%       ^   Prandtl number                  [null]
%
%
%   Examples
%   --------
%       [P,H,S] = refpropsi('PHS','T',300, 'Q', 0.5, 'R134a') ;
%       P = refpropsi('P','T',300, 'Q', 0.5, 'R134a') ;
%       [P,H,S] = refpropsi('PHS','T',300, 'Q', 0.5, {'R134a', 'R152A'}, [0.95, 0.05]) ;
%       props = refpropsi('PHS','T',300, 'Q', 0.5, 'R134a', 1, 'struct') ;
%            -> props.P = Pval, props.H = Hval, props.S = Sval
%       props = refpropsi('PHS','T',300, 'Q', 0.5, 'R134a', 1, 'array') ;
%            -> props = [Pval, Hval, pSval]
%       props = refpropsi('PHS','T',300, 'Q', 0.5, 'R134a', 1, 'cell') ;
%            -> props = {Pval, Hval, pSval}
%       [Pc,Tc] = refpropsi('PT','C',0, ' ', 0, 'R134a') ;
%

    %#ok<*NASGU>
    %#ok<*ASGLU>
    
    %== Variables ==%
    persistent current_fluid ;
    persistent RP_PATH
    RP_CHAR_LEN = 256 ;
    RP_REFS_LEN = 4 ;
    RP_MERR_LEN = 256 ;
    RP_NCOMP_MAX = 20 ;
    LIB_NAME = 'refpropsi' ;

    ierr = 0 ;
    merr = zeros(1, RP_MERR_LEN) ;

    
    %== Check and parse input ==%
    % Check input
    if ~iscell(fluid)
        fluid = {fluid} ;
    end
    if nargin<7
        x = 1 ;
    end
    if nargin<8
        out_type = '' ;
    end
    if nargin<9
        warning_flag = false ;
    end
    if length(x)~=length(fluid)
        error('Refprop:massFraction', 'Fluid and x must have same size') ;
    end
    if abs(sum(x)-1)>1e-10
        error('Refprop:massFraction', 'Sum of x must be 1') ;
    end
    
    % Order input to have always prop name sorted
    prop1 = lower(prop1) ;
    prop2 = lower(prop2) ;
    if prop1<prop2
        prop = [prop1, prop2] ;
        value = [value1; value2] ;
    else
        prop = [prop2, prop1] ;
        value = [value2; value1] ; 
    end
    prop_req_ = prop_req ;
    prop_req = lower(prop_req) ;
    
    
    
    %== Load and setup ==%
    % Find refprop path if needed
    if isempty(RP_PATH)
        RP_PATH = find_refprop_path() ;
    end
    
    % Load refprop dll if needed
    if ~libisloaded(LIB_NAME)
        load_refprop() ;
    end

    % Set fluid if needed
    if ~isequal(current_fluid, fluid)
        current_fluid = fluid ;
        set_fluid(fluid) ;
        check_error()
    end

    
    
    %== Prepare computations ==%
    % Convert mass fraction to molar fraction
    [~, x, ~] = calllib(LIB_NAME, 'XMOLEdll', x, zeros(1,length(x)), 0) ;
    % Compute molar mass
    [~, molm] = calllib(LIB_NAME, 'WMOLdll', x, 0) ;
    molm = molm*1e-3 ;
    
    
    % Critical point if requested
    if strcmp(prop, ' c')
        [x,t,p,d,ierr,merr] = calllib(LIB_NAME, 'CRITPdll', x, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
        check_error() ;
        prop = 'dt' ;
        value = [d*(molm*1e3),t] ;
    end
    
    %== Call flash function ==%
    % Generic flash
    zx = zeros(1, RP_NCOMP_MAX) ;
    if strcmp(prop, 'du')
        [d,e,x,t,p,dl,dv,xl,xv,q,h,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'DEFLSHdll', value(1)/(molm*1e3), value(2)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'dh')
        [d,h,x,t,p,dl,dv,xl,xv,q,e,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'DHFLSHdll', value(1)/(molm*1e3), value(2)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'ds')
        [d,s,x,t,p,dl,dv,xl,xv,q,e,h,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'DSFLSHdll', value(1)/(molm*1e3), value(2)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'su')
        [e,s,x,t,p,d,dl,dv,xl,xv,q,h,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'ESFLSHdll', value(2)*molm, value(1)*molm, x, 0, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'hs')
        [h,s,x,t,p,d,dl,dv,xl,xv,q,e,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'HSFLSHdll', value(1)*molm, value(2)*molm, x, 0, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'dp')
        [p,d,x,t,dl,dv,xl,xv,q,e,h,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'PDFLSHdll', value(2)*1e-3, value(1)/(molm*1e3), x, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'pu')
        [p,e,x,t,d,dl,dv,xl,xv,q,h,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'PEFLSHdll', value(1)*1e-3, value(2)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'hp')
        [p,h,x,t,d,dl,dv,xl,xv,q,e,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'PHFLSHdll', value(2)*1e-3, value(1)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'pq')
        [p,q,x,~,t,d,dl,dv,xl,xv,e,h,s,cv,cp,w,ierr,merr] = calllib(LIB_NAME,'PQFLSHdll', value(1)*1e-3, value(2), x, 2, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'ps')
        [p,s,x,t,d,dl,dv,xl,xv,q,e,h,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'PSFLSHdll', value(1)*1e-3, value(2)*molm, x, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'dt')
        [t,d,x,p,dl,dv,xl,xv,q,e,h,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'TDFLSHdll', value(2), value(1)/(molm*1e3), x, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'tu')
        [t,e,x,~,p,d,dl,dv,xl,xv,q,h,s,cv,cp,w,ierr,merr] = calllib(LIB_NAME,'TEFLSHdll', value(1), value(2)*molm, x, 1, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'ht')
        [t,h,x,~,p,d,dl,dv,xl,xv,q,e,s,cv,cp,w,ierr,merr] = calllib(LIB_NAME,'THFLSHdll', value(2), value(1)*molm, x, 1, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'pt')
        [t,p,x,d,dl,dv,xl,xv,q,e,h,s,cv,cp,w,ierr,merr] =   calllib(LIB_NAME,'TPFLSHdll', value(2), value(1)*1e-3, x, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'qt')
        [t,q,x,~,p,d,dl,dv,xl,xv,e,h,s,cv,cp,w,ierr,merr] = calllib(LIB_NAME,'TQFLSHdll', value(2), value(1), x, 2, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    elseif strcmp(prop, 'st')
        [t,s,x,~,p,d,dl,dv,xl,xv,q,e,h,cv,cp,w,ierr,merr] = calllib(LIB_NAME,'TSFLSHdll', value(2), value(1)*molm, x, 1, 0, 0, 0, 0, zx, zx, 0, 0, 0, 0, 0, 0, ierr, merr, RP_MERR_LEN) ;
    else
        error('Refprop:propertyCombination', 'Unsupported combination: %s', prop) ;
    end
    
    if ierr~=0
        check_error() ;
    end
    
    %== Other properties ==%
    % Transport properties
    if ~isempty(strfind(prop_req, 'v')) || ~isempty(strfind(prop_req, 'l')) || ~isempty(strfind(prop_req, '^'))
        [t,d,x,mu,lambda,ierr,merr] = calllib(LIB_NAME,'TRNPRPdll', t, d, x, 0, 0, ierr, merr, RP_MERR_LEN) ;
        check_error() ;
        if q>0 && q<1
            warning('Refprop:transportProperties', 'Transport properties are not valid with two phases') ;
        end
    end
    
    
    %== Build output ==%
    out = cell(size(prop_req)) ;
    for k = 1:length(prop_req)
        switch prop_req(k)
            case 'a'
                out{k} = w ;
            case 'c'
                out{k} = cp/molm ;
            case 'd'
                out{k} = d*(molm*1e3) ;
            case 'h'
                out{k} = h/molm ;
            case 'i'
                [t,dl,dv,xl,xv,sigma,ierr,merr] = calllib(LIB_NAME, 'SURTENdll', t, dl, dv, xl, xv, 0, ierr, merr, RP_MERR_LEN) ;
                check_error() ;
                out{k} = sigma ;
            case 'k'
                out{k} = cp/cv ;
            case 'l'
                out{k} = lambda ;
            case 'm'
                out{k} = molm ;
            case 'o'
                out{k} = cv/molm ;
            case 'p'
                out{k} = p*1000 ;
            case 'q'
                out{k} = q ;
            case 's'
                out{k} = s/molm ;
            case 't'
                out{k} = t ; 
            case 'u'
                out{k} = e/molm ;
            case 'v'
                out{k} = mu*1e-6 ; 
            case 'z'
                R0 = 8.314459848 ;
                out{k} = p/(d*R0*t) ;
            case '^'
                out{k} = (cp/molm)*(mu*1e-6)/lambda ;
            otherwise
                error('Refprop:property', 'Unsupported property: %s', prop_req(k));
        end
    end
    
    % Change output type if required
    if strcmp(out_type, 'array')
        varargout = {cell2mat(out)} ;
    elseif strcmp(out_type, 'cell')
        varargout = {out} ;
    elseif strcmp(out_type, 'struct')
        varargout = struct() ;
        for k = 1:length(out)
            if prop_req_(k)=='^'
                varargout.Pr = out{k} ;
            else
                varargout.(prop_req_(k)) = out{k} ;
            end
        end
        varargout = {varargout} ;
    else
        varargout = out ;
    end
    
    %=====================================================================%
    function load_refprop()
        % Load refprop dll
        arch = computer('arch') ;

        % Check computer architecture
        switch(arch)
            case 'win32'
                dllname = 'refprop.dll' ;
                prototype = @() rp_proto(RP_PATH) ;
                
            case 'win64'
                dllname = 'refprp64.dll' ;
                prototype = @() rp_proto64(RP_PATH) ;
                
            case 'glnxa64'
                dllname = 'librefprop.so'  ;
                prototype = @() rp_proto64(RP_PATH) ;
                
            otherwise
                error('Refprop:unsupported', 'Unsupported architecture: %s', arch) ;
        end

        % Load refprop
        path = fullfile(RP_PATH, dllname) ;
        loadlibrary(path, prototype, 'alias', LIB_NAME) ;
    end

    function set_fluid(fluid)
        % Set fluid in refprop
        nf = length(fluid) ;
        
        fluidstr = cellfun(@(s) fullfile(RP_PATH, 'fluids', [s,'.fld']), fluid, 'UniformOutput', 0) ;
        fluidstr = strjoin(fluidstr,'|');
        fluidstr = [fluidstr, zeros(1, RP_CHAR_LEN*RP_NCOMP_MAX-length(fluidstr))] ;
        
        mixstr = fullfile(RP_PATH, 'fluids', 'hmx.bnc');
        mixstr = [mixstr, zeros(1, RP_CHAR_LEN*RP_NCOMP_MAX-length(mixstr))] ;
        
        refstr = ['DEF', 0] ;
        
        [~,~,~,~,ierr,merr] = calllib(LIB_NAME, 'SETUPdll', nf, int32(fluidstr), int32(mixstr), int32(refstr), ierr, merr, RP_CHAR_LEN*RP_NCOMP_MAX, RP_CHAR_LEN, RP_REFS_LEN, RP_MERR_LEN) ;
    end

    function path = find_refprop_path()
        % Find refprop folder
        arch = computer('arch') ;

        % Check computer architecture
        path = '' ;
        switch(arch)
            case {'win32', 'win64'}
                path = 'C:\Program files\Refprop\' ;
                if ~exist(path, 'dir')
                    path = 'C:\Program files (x86)\Refprop\' ;
                end
                
            case 'glnxa64'
                path = '/opt/refprop/' ;
                
            otherwise
                error('Refprop:unsupported', 'Unsupported architecture: %s', arch) ;
        end
        
        if ~exist(path, 'dir')
            error('Refprop:path', 'Unable to find refprop folder') ;
        end
    end

    function check_error()
        % Check if an error occured in Refprop
        if ierr>0
            error('Refprop:dllError', char(merr)) ;
        elseif warning_flag && ierr<0
            warning('Refprop:dllWarning', char(merr)) ;
        end
    end
end
