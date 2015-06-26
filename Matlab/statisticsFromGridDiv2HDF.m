function res = statisticsFromGridDiv2HDF(basename,first,last,lBeg,lEnd,threshold,run);

fieldname_list = 'lambda1,lambda2,lambda3,sss,wws,s2,u2,w2,cosss,cosww,coswl1,coswl2,coswl3,divref,ii,cosulambda1m,cosulambda2m,cosulambda3m';
dfields = {'lambda1','lambda2','lambda3','sss','wws','s2','u2','w2','cosss','cosww','coswl1','coswl2','coswl3','divref','ii','cosulambda1m','cosulambda2m','cosulambda3m'};
dtype = repmat({'double'},18,1);


if ~nargin,
    basename = {'F:/PTV/convection/150/res_a/grid.',...
        'F:/PTV/convection/150/res_b/grid.',...
        'F:/PTV/convection/150/res_c/grid.',...
        'F:/PTV/convection/150/res_d/grid.',...
        'F:/PTV/convection/150/res_e/grid.',...
        'F:/PTV/convection/150/res_f/grid.';
    'F:/convection/newRun_90/res_2707a/grid.',...
        'F:/convection/newRun_90/res_2707b/grid.',...
        'F:/convection/newRun_90/res_2707c/grid.',...
        'F:/convection/newRun_90/res_2707d/grid.',...
        'F:/convection/newRun_90/res_2707e/grid.',...
        'F:/convection/newRun_90/res_2707f/grid.'};
   
    first = 10200;
    last = 10300;
    lBeg = 1;
    lEnd = 1;
    threshold = 1e-6;
    offset = 10;
    run = 1;
end



% CounterN=1;
% CounterP=1;
stat.wws=zeros(10000,1);
stat.sss=zeros(10000,1);
stat.coswWArray=zeros(10000,1);
stat.coswl1=zeros(10000,1);
stat.coswl2=zeros(10000,1);
stat.coswl3=zeros(10000,1);
stat.ei=zeros(10000,3);
stat.Ssq=zeros(10000,1);
stat.Energy=zeros(10000,1);
stat.cosS=zeros(10000,1);
stat.Wsq=zeros(10000,1);
stat.divref=zeros(10000,1);
stat.ii=zeros(10000,1);
stat.cosulambdam = zeros(10000,3);



% load meanField;
% load meanCorField;

if lBeg == lEnd
    load(['meanField_',int2str(run),'_',int2str(lBeg)]);
    load(['meanCorField_',int2str(run),'_',int2str(lBeg)]);
    hfile = ['meanCorField_',int2str(run),'_',int2str(lBeg),'.hdf']
else
    load(['meanField'])
    load(['meanCorField'])
    hfile = ['meanCorField','.hdf']
end

try 
    if hdf('H','ishdf', hfile)
        access = 'write';
        %     answer = input('Do you want to overwrite ... (y/n)  ','s');
        %     if strcmp(answer,'y')
        %       hdfml('closeall'), delete(hfile), access = 'create';
        %       feval('TrajPoint2HDF_v7',hfile,varargin{:});
        %     else
        warning('File exists, delete it manually') % this option should be changed, if somebody wants to append 28.11.03
        return
        % end
    else
        access = 'create';
        % open new or existing hfile
        file_id = hdf('H','open', hfile, access, 0);
        if file_id == -1
            error('HDF hopen failed');
        end
        
        % initialize the V interface
        status = hdf('V','start',file_id);
        if status == -1
            error('HDF vstart failed');
        end
    end
    
    runIndx = 0;
    
    for loop = lBeg:lEnd
        for i = first:last
            runIndx = runIndx + 1;
            clear f;
            [run loop i]
            f = load([basename{run,loop},int2str(i)]);
            s = size(f);
            if s(1,1) > 0
                u = f(:,4)-mea.u;
                v = f(:,5)-mea.v;
                w = f(:,6)-mea.w;
                ux = f(:,7)-mea.ux;
                uy = f(:,8)-mea.uy;
                uz = f(:,9)-mea.uz;
                vx = f(:,10)-mea.vx;
                vy = f(:,11)-mea.vy;
                vz = f(:,12)-mea.vz;
                wx = f(:,13)-mea.wx;
                wy = f(:,14)-mea.wy;
                wz = f(:,15)-mea.wz;
                w1 = wy-vz;
                w2 = uz-wx;
                w3 = vx-uy;
                absw = (w1.^2+w2.^2+w3.^2).^0.5;
                enstrophy = w1.^2+w2.^2+w3.^2;
                s11 = 0.5*(ux+ux);
                s12 = 0.5*(uy+vx);
                s13 = 0.5*(uz+wx);
                s22 = 0.5*(vy+vy);
                s23 = 0.5*(vz+wy);
                s33 = 0.5*(wz+wz);
                
                wS1 = w1.*s11+w2.*s12+w3.*s13;
                wS2 = w1.*s12+w2.*s22+w3.*s23;
                wS3 = w1.*s13+w2.*s23+w3.*s33;
                absws = (wS1.^2+wS2.^2+wS3.^2).^0.5;
                
                divref = f(:,31);
                
                
                % ------------- HDF part -------------------
                % create a new vs data 
                access  =  'w';
                vdata_ref  =  -1;                                         % flag to create
                vdata_id  =  hdf('VS','attach', file_id, vdata_ref, access);
                %         vdata_id  =  hdf('V',s('attach', file_id, vdata_ref, access);
                if vdata_id  ==  -1
                    error('HDF vdata attach failed')
                end
                
                %  offset  =  0;
                
                
                % assign the group the name "i"
                % sname  =  sprintf('%d', lastTrajIndx);
                sname  =  sprintf('%d', runIndx);
                status  =  hdf('VS','setname', vdata_id, sname);
                if status  ==  -1
                    error('HDF vsetname failed')
                end
                
                % assign the group the class "Trajectory"
                sclass  =  'Grid';
                status  =  hdf('VS','setclass', vdata_id, sclass);
                if status  ==  -1
                    error('HDF vsetclass failed')
                end
                
                for j  =  1 : length(dfields)                 % 32 fields
                    status  =  hdf('VS','fdefine',vdata_id, dfields{j},'double',1);
                end
                
                status  =  hdf('VS','setfields',vdata_id,fieldname_list);
                if status  ==  -1
                    error('HDF VS setfields failed')
                end
                
                % Interlace mode is set to full, meaning that data is written row
                % by row of the record, rather than by field (not all x, then all
                % y, but all x,y,z, ... for the record 1, then for 2, etc)
                status  =  hdf('VS','setinterlace',vdata_id,'full');
                if status  ==  -1
                    error('HDF VS setinterlace failed')
                end
                
                
                CounterN = 0;
                
                for ii = 1:s(1,1)
                    A = [s11(ii) s12(ii) s13(ii); s12(ii) s22(ii) s23(ii); s13(ii) s23(ii) s33(ii)];
                    if divref(ii) < 0.1 %& all(isfinite(A(:))) 
                       %  A = [s11(ii) s12(ii) s13(ii); s12(ii) s22(ii) s23(ii); s13(ii) s23(ii) s33(ii)];
                        
                        [V,D] = eig(A);
                        [D,k]  =  sort(diag(D)); 	% ascending order
                        D  =  diag(D(end:-1:1)); 			% descending order		
                        V  =  V(:,k(end:-1:1)); 	% the same order for eigenvectors
                        
                        %                         if D(1,1)^2 + D(2,2)^2 + D(3,3)^2 < 100
                        
                        CounterN = CounterN+1
                        
                        stat.ei(CounterN,1) = D(1,1);
                        stat.ei(CounterN,2) = D(2,2);
                        stat.ei(CounterN,3) = D(3,3);
                        
                        
                        stat.sss(CounterN,1)  =  -D(1,1)*D(2,2)*D(3,3);                        
                        stat.wws(CounterN,1)  =  ([w1(ii) w2(ii) w3(ii)]*[wS1(ii) wS2(ii) wS3(ii)]');
                        %                             stat.Ssq(CounterN,1)  =  mea.D1(ii)^2+mea.D2(ii)^2+mea.D3(ii)^2;
                        stat.Ssq(CounterN,1)  =  D(1,1)^2+D(2,2)^2+D(3,3)^2;
                        %                             stat.Energy(CounterN,1) = 0.5*(mea.u(ii)^2+mea.v(ii)^2+mea.w(ii)^2);
                        stat.Energy(CounterN,1) = 0.5*(u(ii)^2+v(ii)^2+w(ii)^2);
                        stat.Wsq(CounterN,1) = absw(ii)^2;
                        %%%%%%% added by michele Negative%%%%%%%%%%%%%%                        
                        S11 = -s11(ii)^2-s12(ii)^2-s13(ii)^2;
                        S12 = -s11(ii)*s12(ii)-s12(ii)*s22(ii)-s13(ii)*s23(ii);
                        S13 = -s11(ii)*s13(ii)-s12(ii)*s23(ii)-s13(ii)*s33(ii);
                        S23 = -s12(ii)*s13(ii)-s22(ii)*s23(ii)-s23(ii)*s33(ii);
                        S22 = -s22(ii)^2-s12(ii)^2-s23(ii)^2;
                        S33 = -s33(ii)^2-s13(ii)^2-s23(ii)^2;                                               
                        S2 = S11.^2+S22.^2+S33.^2+2*(S12.^2+S13.^2+S23.^2);                        
                        cosSs = -D(1,1)*D(2,2)*D(3,3)/((D(1,1)^2+D(2,2)^2+D(3,3)^2)^0.5 * S2^0.5);                        
                        stat.cosS(CounterN,1) = cosSs;
                        
                        
                        %%%%%%% added by michele Negative%%%%%%%%%%%%%%    
                          
                         lam1m=[mea.V1x(ii) mea.V1y(ii) mea.V1z(ii)];
                         lam2m=[mea.V2x(ii) mea.V2y(ii) mea.V2z(ii)];
                         lam3m=[mea.V3x(ii) mea.V3y(ii) mea.V3z(ii)];
                        
                         magn1=(lam1m(1)^2+lam1m(2)^2+lam1m(3)^2)^0.5;
                         magn2=(lam2m(1)^2+lam2m(2)^2+lam2m(3)^2)^0.5;
                         magn3=(lam3m(1)^2+lam3m(2)^2+lam3m(3)^2)^0.5;
                         
                         magnu=(u(ii)^2+v(ii)^2+w(ii)^2)^0.5;
                         
                         cosuLam1=([u(ii) v(ii) w(ii)]*[lam1m(1) lam1m(2) lam1m(3)]')/(magn1*magnu);
                         cosuLam2=([u(ii) v(ii) w(ii)]*[lam2m(1) lam2m(2) lam2m(3)]')/(magn2*magnu);
                         cosuLam3=([u(ii) v(ii) w(ii)]*[lam3m(1) lam3m(2) lam3m(3)]')/(magn3*magnu);
                        
                        stat.cosulambdam(CounterN,1)=cosuLam1;
                        stat.cosulambdam(CounterN,2)=cosuLam2;
                        stat.cosulambdam(CounterN,3)=cosuLam3;
                        
                        
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        coswW = ([w1(ii) w2(ii) w3(ii)]*[wS1(ii) wS2(ii) wS3(ii)]')/(absw(ii)*absws(ii));
                        hiwl1 = ([w1(ii) w2(ii) w3(ii)]*V(:,1))/absw(ii);
                        hiwl2 = ([w1(ii) w2(ii) w3(ii)]*V(:,2))/absw(ii);
                        hiwl3 = ([w1(ii) w2(ii) w3(ii)]*V(:,3))/absw(ii);
                        stat.coswWArray(CounterN,1) = coswW;
                        stat.coswl1(CounterN,1) = hiwl1;
                        stat.coswl2(CounterN,1) = hiwl2;
                        stat.coswl3(CounterN,1) = hiwl3;
                        
                        stat.divref(CounterN,1) = divref(ii);
                        stat.ii(CounterN,1) = ii;
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        %                         end
                    end
                end
                
                
                
                
                % write the buffered data into the first vdata with full interlace mode
                num_of_records  =  hdf('VS','write',vdata_id,{stat.ei(1:CounterN,1)',stat.ei(1:CounterN,2)',stat.ei(1:CounterN,3)',...
                        stat.sss(1:CounterN,1)',stat.wws(1:CounterN,1)',stat.Ssq(1:CounterN,1)',stat.Energy(1:CounterN,1)',...
                        stat.Wsq(1:CounterN,1)',stat.cosS(1:CounterN,1)',stat.coswWArray(1:CounterN,1)',stat.coswl1(1:CounterN,1)',...
                        stat.coswl2(1:CounterN,1)',stat.coswl3(1:CounterN,1)',stat.divref(1:CounterN,1)',stat.ii(1:CounterN,1)',...
                    stat.cosulambdam(1:CounterN,1)',stat.cosulambdam(1:CounterN,2)',stat.cosulambdam(1:CounterN,3)'});
                
               
                
                % detach from the first vdata
                status  =  hdf('VS','detach',vdata_id);
                if status  ==  -1, error('HDF vsdetach failed'), end
                %         lastTrajIndx  =  lastTrajIndx + 1; 
                %     end  % of all trajectories in the file
                
            end % for s(1,1)
        end % first:last
    end % loop
    % ---------------------------------------------------------------------------------------------------------        
    
    
    % end vgroup interface access
    status  =  hdf('V','end',file_id);
    if status  ==  -1
        error('HDF v end failed')
    end
    %     
    %     % end vgroup interface access
    %     status  =  hdf('SD','end',sd_id); 
    %     if status  ==  -1
    %         error('HDF sd end failed')
    %     end
    
    % close the HDF file
    status  =  hdf('H','close',file_id);
    if status  ==  -1
        error('HDF hclose failed')
    end
    % ---------------------------------------------------------------------------------------------------------
    
    
    %     cd(wd);
    
catch
    hdfml('closeall')
    %     cd(wd)
    lasterr
end


