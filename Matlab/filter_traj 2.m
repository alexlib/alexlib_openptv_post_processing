function ptv = filter_traj(traj, ptv)

% global traj;
% global ptv;

for i=1:length(traj)
    if i==1
        nx(i)=0.5*ptv(traj(i,1),traj(i,2),3)+0.5*ptv(traj(i+1,1),traj(i+1,2),3);
        ny(i)=0.5*ptv(traj(i,1),traj(i,2),4)+0.5*ptv(traj(i+1,1),traj(i+1,2),4);
        nz(i)=0.5*ptv(traj(i,1),traj(i,2),5)+0.5*ptv(traj(i+1,1),traj(i+1,2),5);
        
        
        ntot_pix(i)=0.5*ptv(traj(i,1),traj(i,2),6)+0.5*ptv(traj(i+1,1),traj(i+1,2),6);
        nx_pix(i)=0.5*ptv(traj(i,1),traj(i,2),7)+0.5*ptv(traj(i+1,1),traj(i+1,2),7);
        ny_pix(i)=0.5*ptv(traj(i,1),traj(i,2),8)+0.5*ptv(traj(i+1,1),traj(i+1,2),8);
        ngrv(i)=0.5*ptv(traj(i,1),traj(i,2),9)+0.5*ptv(traj(i+1,1),traj(i+1,2),9);
        
    end
    if i==length(traj)
        nx(i)=0.5*ptv(traj(i,1),traj(i,2),3)+0.5*ptv(traj(i-1,1),traj(i-1,2),3);
        ny(i)=0.5*ptv(traj(i,1),traj(i,2),4)+0.5*ptv(traj(i-1,1),traj(i-1,2),4);
        nz(i)=0.5*ptv(traj(i,1),traj(i,2),5)+0.5*ptv(traj(i-1,1),traj(i-1,2),5);
        
        
        ntot_pix(i)=0.5*ptv(traj(i,1),traj(i,2),6)+0.5*ptv(traj(i-1,1),traj(i-1,2),6);
        nx_pix(i)=0.5*ptv(traj(i,1),traj(i,2),7)+0.5*ptv(traj(i-1,1),traj(i-1,2),7);
        ny_pix(i)=0.5*ptv(traj(i,1),traj(i,2),8)+0.5*ptv(traj(i-1,1),traj(i-1,2),8);
        ngrv(i)=0.5*ptv(traj(i,1),traj(i,2),9)+0.5*ptv(traj(i-1,1),traj(i-1,2),9);
    end
    if i>1 & i<length(traj)
        order=min(i-1,length(traj)-i);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if order>4
            order=4;
        end
        switch order
            case 1
                weight=[2 1]; % Here weight is defined by Pascal's triangle
            case 2
                weight=[6 4 1];
            case 3
                weight=[20 15 6 1];
            case 4
                weight=[70 56 28 8 1];
                
            case 5
                weight=[252 210 120 45 10 1];
                
            case 6
                weight=[924 792 495 220 66 12 1];
                
            case 7
                weight=[3432,3003,2002,1001,364,91,14,1];
                
            case 8
                weight=[12870,11440,8008,4368,1820,560,120,16,1];
            case 9
                weight=[48620,43758,31824,18564,8568,3060,816,153,18,1];
            case 10
                weight=[184756,167960,125970,77520,38760,15504,4845,1140,190,20,1];
            case 11
                weight=[705432,646646,497420,319770,170544,74613,26334,7315,1540,231,22,1];
            case 12
                weight=[2704156,2496144,1961256,1307504,735471,346104,134596,42504,10626,2024,276,24,1];
                
            case 13
                weight=[10400600,9657700,7726160,5311735,3124550,1562275,657800,230230,65780,14950,2600,325,26,1];
                
            case 14
                weight=[40116600,37442160,30421755,21474180,13123110,6906900,3108105,1184040,376740,98280,20475,3276,378,28,1];
                
            case 15
                weight=[155117520,145422675,119759850,86493225,54627300,30045015,14307150,5852925,2035800,593775,142506,27405,4060,435,30,1];
                
            case 16
                weight=[601080390,565722720,471435600,347373600,225792840,129024480,64512240,28048800,10518300,3365856,906192,201376,35960,4960,496,32,1];
                
            case 17
                weight=[2333606220.00000,2203961430.00000,1855967520.00000,1391975640.00000,927983760,548354040,286097760,131128140,52451256,18156204,5379616,1344904,278256,46376,5984,561,34,1];
            case 18
                weight=[9075135300.00000,8597496600.00000,7307872110.00000,5567902560.00000,3796297200.00000,2310789600.00000,1251677700.00000,600805296,254186856,94143280,30260340,8347680,1947792,376992,58905,7140,630,36,1];
            case 19
                weight=[35345263800.0000,33578000610.0000,28781143380.0000,22239974430.0000,15471286560.0000,9669554100.00000,5414950296.00000,2707475148.00000,1203322288.00000,472733756,163011640,48903492,12620256,2760681,501942,73815,8436,703,38,1];
            case 20
                weight=[137846528820.000,131282408400.000,113380261800.000,88732378800.0000,62852101650.0000,40225345056.0000,23206929840.0000,12033222880.0000,5586853480.00000,2311801440.00000,847660528,273438880,76904685,18643560,3838380,658008,91390,9880,780,40,1];
        end
        su=0;
        fx=0;fy=0;fz=0;ftotpix=0;fxpix=0;fypix=0;fgrv=0;
        for j=i-order:i+order
            su=su+weight(abs(i-j)+1);
            fx=fx+ptv(traj(j,1),traj(j,2),3)*weight(abs(i-j)+1);
            fy=fy+ptv(traj(j,1),traj(j,2),4)*weight(abs(i-j)+1);
            fz=fz+ptv(traj(j,1),traj(j,2),5)*weight(abs(i-j)+1);
            
            ftotpix=ftotpix+ptv(traj(j,1),traj(j,2),6)*weight(abs(i-j)+1);
            fxpix=fxpix+ptv(traj(j,1),traj(j,2),7)*weight(abs(i-j)+1);
            fypix=fypix+ptv(traj(j,1),traj(j,2),8)*weight(abs(i-j)+1);
            fgrv=fgrv+ptv(traj(j,1),traj(j,2),9)*weight(abs(i-j)+1);
            
            
        end
        nx(i)=fx/su;
        ny(i)=fy/su;
        nz(i)=fz/su;
        
        ntot_pix(i)=ftotpix/su;
        nx_pix(i)=fxpix/su;
        ny_pix(i)=fypix/su;
        ngrv(i)=fgrv/su;
        
        
        
        
    end
end

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),3)=nx(i);
    ptv(traj(i,1),traj(i,2),4)=ny(i);
    ptv(traj(i,1),traj(i,2),5)=nz(i);
    
    ptv(traj(i,1),traj(i,2),6)=ntot_pix(i);
    ptv(traj(i,1),traj(i,2),7)=nx_pix(i);
    ptv(traj(i,1),traj(i,2),8)=ny_pix(i);
    ptv(traj(i,1),traj(i,2),9)=ngrv(i);
    
    
end
