function [SLD, DATA, REF] = SLD_function(DATA, REF, pars)
% function [SLD] = SLD_function(DATA, REF, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : SLD technique withouth SNR evaluation
% INPUTS: 
%         - DATA : Sample data
%                 - DATA.Bmode
%                 - DATA.rf
%                 - DATA.fc
%                 - DATA.fs
%                 - DATA.x
%                 - pars.z
%         - REF  : Reference phantom
%                 - REF.Bmode
%                 - REF.rf
%                 - REF.fc
%                 - REF.fs
%                 - REF.x
%                 - REF.z
%         - pars : 
%                 - pars.cs % speed of sound SoS in [m/s]
%                 - pars.z_ini; % ROI ini axial (depth) % in [m]
%                 - pars.z_end; % ROI end axial (depth) % in [m]
%                 - pars.x_ini; % ROI ini lateral % in [m]
%                 - pars.x_end; % ROI ini lateral % in [m]
%                 - pars.bw BANDWITH i.e [ 4 10 ] % in [MHz]
%                 - pars.nb_lambda_axial % datablock in wavelengths
%                 - pars.overlap_axial; 
%                 - pars.overlap_lateral;
%                 - pars.REF_num % CHOICE OF REFERENCE PHANTOM CHECK attenuation_phantoms_Np.m 
%                 - pars.SLOPE % In case of linear attenuation (simulation)
% OUTPUTS: 
%         - SLD : 
%                 - SLD_term: Spectral Log Difference Term of nºchannels of frequency 3D array (axial, lateral, band)
%                 - band: Frequency band vector [MHz]
%                 - A : A matrix por SLD as 2D-Inverse problem [-4 zp-zd f, speye(m*n)  | 1*f speye(m*n)]
%                 - x_ori : DATA.x 
%                 - z_ori : DATA.z 
%                 - x : X_ROI
%                 - z : Z_ROI
%                 - SLD.Bmode = DATA.Bmode;
% AUTHORs: Edmundo Arom Miranda & Andres Coila, based on LIM repository 
%          (Andres Coila Code RSLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                saran_layer     = 0;

                sam1 = DATA.rf(:,:);
                rf1 = REF.rf1(:,:);
                rf2 = REF.rf2(:,:);
                rf3 = REF.rf3(:,:);
                rf4 = REF.rf4(:,:);

                fs = REF.fs;
                x = REF.x;
                z = REF.z;

                z_inf = pars.z_ini; % in [m]
                z_sup = pars.z_end; % in [m]
                x_inf = pars.x_ini;  % in [m];
                x_sup = pars.x_end; % in [m];


                freq_L = pars.bw(1);
                freq_H = pars.bw(2);
                blocksize = pars.nb_lambda_axial;
                overlap_pc = pars.overlap_axial;
                overlap_pc = pars.overlap_lateral;
                window_type = pars.window_type;
                c0 = pars.cs;
                REF_num = pars.REF_num;

                Imin  = find(DATA.z > z_inf);  
                Imax  = find(DATA.z > z_sup);  
            
                Jmin = find(DATA.x > x_inf);  
                Jmax = find(DATA.x > x_sup);

                env_rfdata_sam = abs(hilbert( DATA.rf ));
                env_rfdata_sam_roi = env_rfdata_sam (Imin(1):Imax(1), Jmin(1):Jmax(1), 1);


                
                % Lateral (x) and axial (z) deltas                
                dx=(x(end)-x(1))/(length(x)-1);
                dz=(z(end)-z(1))/(length(z)-1);
                
                % x and z at the starting time
                x_ori = x;
                z_ori = z;
                
                % x and z in [cm]
                x = x;  
                z = z;  
                
                % Plot B-mode image
                font = 14;
                dyn_range = 60;
  
                Bmode = db(hilbert(DATA.rf));
                Bmode = Bmode - max(Bmode(:));

                Im_db_ori = Bmode;
                % Im_db_ori = DATA.Bmode;

                
                if REF_num == 111 % CASE LINEAR DEPENDENCY
                    SLOPE = pars.SLOPE; % [dB/cm-Hz]
                else
                    SLOPE = NaN; % DISCARD
                end

                % Cut until the limits for ACS estimation
                ind_x = x_inf <= x & x <= x_sup;
                x = x(ind_x);
                ind_z = z_inf <= z & z <= z_sup;
                z = z(ind_z);
                sam1 = sam1(ind_z, ind_x);
                
                ref1 = rf1;
                ref2 = rf2;
                ref3 = rf3;
                ref4 = rf4;

                ref1 = ref1(ind_z,ind_x);
                ref2 = ref2(ind_z,ind_x);
                ref3 = ref3(ind_z,ind_x);
                ref4 = ref4(ind_z,ind_x);

                
                freq_L = freq_L*1e6;   % (Hz)
                freq_H = freq_H*1e6;   % (Hz)
                
                wl = c0/mean([freq_L freq_H]);   % Wavelength (m)
                disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
                    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
                
                % Number of repetitions of datablocks within a single
                % datablock
                rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
                rpt = round(rpt);
                overlap = 1 - (1/rpt);
                
                % Part without lateral overlap
                wx  = round((blocksize*wl)/(dx*rpt));   
                % Number of lines laterally = repetitions * block without repetition
                nx  = rpt*wx;   
                new_blocksize = round(nx*dx/(wl));
                
                % rf data columns
                L2   = size(sam1,2);   
                % Blocksize colums
                ncol = floor((L2-(rpt-1)*wx)/wx);     
                sam1 = sam1(:,1:wx*(ncol+rpt-1));
                % Actual rf data columns
                L2 = size(sam1,2);   
                x  = x(1:L2);
                
                xi = 1;
                xf = L2;
                x0 = (xi:wx:xf+1-nx);
                x_ACS = x(1,x0+round(nx/2));
                nx_ACS = nx*dx;
                
                n  = length(x0);
                
                wz = floor(nx*dx/(dz*rpt));
                nz = rpt*wz;
                
                % winsize: Percentage of window (max 0.5)
                % nw: Samples of each window axially
                winsize = 0.5;
                nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;  
                L = (nz - nw)*dz*100;   % (cm)
                                
                NFFT = 2^(nextpow2(nw)+2);
                band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
                rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
                f  = band(rang)*1e-6; % [MHz]
                L3 = length(rang);
                p = L3;
                
                L1   = size(sam1,1);   % rf data: rows
                nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
                sam1 = sam1(1:wz*(nrow+rpt-1),:);
                L1   = size(sam1,1);   % rf data: rows
                z    = z(1:L1);
                
                zi = 1;
                zf = L1;
                z0 = (zi:wz:zf+1-nz);
                m  = length(z0);
                
                z_ACS = z(z0+round(nz/2));
                nz_ACS = nz*dz;
                
                z0p = z0 + (nw-1)/2;
                z0d = z0 + (nz-1) - (nw-1)/2;
                                
                %% Plot region of interest B-mode image

                Im=abs(hilbert(sam1));   % envelope calculation
                Im_db=20*log10(Im/max(Im(:)));   % log scale


               
                    ref1 = ref1(1:L1,1:L2);
                    ref2 = ref2(1:L1,1:L2);
                    ref3 = ref3(1:L1,1:L2);
                    ref4 = ref4(1:L1,1:L2);

                disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
                    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
                disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
                disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
                disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
                disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);
                
                
                %bottom = input('Bottom ACS (ejm: 0): ');
                %top = input('Top ACS (ejm: 1.5): ');
                
                %window_type = input('Windowing type (Tukey-0.25 (5), Hamming (6), Boxcar (7)): ');
%                 switch window_type
%                     case 5
%                         windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
%                     case 6
%                         windowing = hamming(nw);   % Hamming
%                     case 7
%                         windowing = rectwin(nw);   % Boxcar
%                     
%                 end
%                 
%                 % Windowing neccesary before Fourier transform
%                 windowing = windowing*ones(1,nx);   
                windowing = window_choice(length(nx), window_type);
                
                figure(103); set(103,'units','normalized','outerposition',[0 0 1 1]);
                figure(106); set(106,'units','normalized','outerposition',[0 0 1 1]);
                
                for jj=1:n
                    for ii=1:m
                        
                        xw = x0(jj) ;   % x window
                        zp = z0p(ii);
                        zd = z0d(ii);
                        
                        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

                        % aux EMZ for SNR % MODIFICATION 02/01/24,
                        % env_block could be pre-saved out of loop
                        block_dp = sam1(zp-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
                        env_block = abs(hilbert( block_dp ));
                        snr_block = mean(env_block)/std(env_block);
                        og_SNR(ii,jj) = snr_block;
                        delta_SNR(ii,jj) = abs(snr_block - 1.91).*100./1.91;                        
                        % aux EMZ

                        blockP(ii,jj) = effective_lines( sam1(zp-(nw-1)/2:2:zp+(nw-1)/2,xw:xw+nx-1) );
                        blockD(ii,jj) = effective_lines( sam1(zd-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );
                        blockT(ii,jj) = effective_lines( sam1(zp-(nw-1)/2:2:zd+(nw-1)/2,xw:xw+nx-1) );
                        
                        [tempSp,temp_psnr_Sp] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                        [tempSd,temp_psnr_Sd] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
                        
                        Sp(ii,jj,:) = (tempSp(rang));
                        Sd(ii,jj,:) = (tempSd(rang));
                        
                        psnr_Sp(ii,jj) = temp_psnr_Sp;
                        psnr_Sd(ii,jj) = temp_psnr_Sd;
                        
                        if jj==floor(3*n/6)
                            if ii==ceil(m/6)
                                figure(103)
                                set(gca,'FontSize',font);
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k');
                                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font);
                                
                                figure(106)
                                spmax = max(tempSp((1:NFFT/2+1)));
                                set(gca,'FontSize',font);                                
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'k');
                                %title('Spectrum'); 
                                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font); 
                                
%                                 figure(7); set(7,'units','normalized','outerposition',[0 0 1 1]);
%                                 s1 = sub_block_p(:,round(end/2));
%                                 smax = max(abs(s1));
%                                 set(gca,'FontSize',font);
%                                 plot(s1/smax,'k');
%                                 %title('Gated ultrasound signal'); 
%                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
%                                 %axis([0 25 -70 0]);
%                                 set(gca,'FontSize',font);  
%                                 axis tight;
%                                 ylim([-1 1])
                            elseif ii==round(m/2)
                                figure(103)
                                hold on;
                                set(gca,'FontSize',font);
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r');
                                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font);

                                figure(106)
                                hold on;
                                set(gca,'FontSize',font);
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'r');
                                %title('Spectrum');
                                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font);
                                
%                                 figure(8); set(8,'units','normalized','outerposition',[0 0 1 1]);
%                                 s2 = sub_block_p(:,round(end/2));
%                                 %hold on;
%                                 set(gca,'FontSize',font);
%                                 plot(s2/smax,'r');
%                                 %title('Gated ultrasound signal'); 
%                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
%                                 %axis([0 25 -70 0]);
%                                 set(gca,'FontSize',font);  
%                                 axis tight;
%                                 ylim([-1 1])
                                
                            elseif ii==floor(4.5*m/6)
                                figure(103)
                                hold on;
                                set(gca,'FontSize',font);
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b');   %
                                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font), 
                                grid on;
                                %pause
                                legend('Top','Half','Bottom');

                                figure(106)
                                hold on;
                                set(gca,'FontSize',font);
                                plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/spmax),'b');   %
                                %title('Spectrum'); 
                                xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity normalized (dB)');
                                axis([0 25 -70 0]);
                                set(gca,'FontSize',font); grid minor;
                                %pause
                                legend('1 cm above focal dist.','At focal distance','1 cm below focal dist.');
                                
%                                 figure(9); set(9,'units','normalized','outerposition',[0 0 1 1]);
%                                 %hold on;
%                                 s3 = sub_block_p(:,round(end/2));
%                                 set(gca,'FontSize',font);
%                                 plot(s3/smax,'b');
%                                 %title('Gated ultrasound signal'); 
%                                 xlabel('\bfSamples'); ylabel('\bfVoltage normalized (V)');
%                                 %axis([0 25 -70 0]);
%                                 set(gca,'FontSize',font);  
%                                 %legend('Top','Half','Bottom');
%                                 axis tight;
%                                 ylim([-1 1])                    
                                
                            end
                            
                        end
                        
                    end
                    
                end
                figure(106)
                hold on;                
                plot(band((1:NFFT/2+1))*1e-6,-20*ones(NFFT/2+1,1),'k--');
                plot(band((1:NFFT/2+1))*1e-6,-15*ones(NFFT/2+1,1),'k--');
                plot(band((1:NFFT/2+1))*1e-6,-10*ones(NFFT/2+1,1),'k--');

                figure(103)
                hold on;                
                title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                       
                plot(band((1:NFFT/2+1))*1e-6,-20*ones(NFFT/2+1,1),'k--');
                plot(band((1:NFFT/2+1))*1e-6,-15*ones(NFFT/2+1,1),'k--');
                plot(band((1:NFFT/2+1))*1e-6,-10*ones(NFFT/2+1,1),'k--');
                
                %% Correlation in reference phantom
                blockP_avg = mean2(blockP);
                blockD_avg = mean2(blockD);
                %blockT_avg = mean2(blockT)
                rxx = 5;
                %mxx = floor((1-overlap)*size(blockT,1));
                %nxx = floor((1-overlap)*size(blockT,2));
                index_block_row = (1:rxx:size(blockT,1));
                index_block_col = (1:rxx:size(blockT,2));
                %blockT2 = blockT(index_block_row,index_block_col);
                %blockT2_avg = mean2(blockT2)
                        
                att_ref = attenuation_phantoms_Np(f, REF_num, SLOPE); % [Np/cm]
                
                %% Diffraction compensation
      
                % Attenuation reference
%                 for jj=1:n
%                     for ii=1:m
%                         att_ref_map(ii,jj,:) = att_ref;
%                     end
%                 end     

                % OPTMIZATION WAY
                att_ref_map = ones(m, n, length(att_ref)) .* reshape(att_ref, 1, 1, []);
                        
                        % Four 'samples' of the reference phantom
                        
                        for jj=1:n
                            for ii=1:m                                
                                xw = x0(jj) ;   % x window
                                zp = z0p(ii);
                                zd = z0d(ii);                                
                                % Reference 1
                                sub_block_p = ref1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                                sub_block_d = ref1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
                                [tempSp1,temp_psnr_Sp1] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                                [tempSd1,temp_psnr_Sd1] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
                                % Reference 2
                                sub_block_p = ref2(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                                sub_block_d = ref2(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
                                [tempSp2,temp_psnr_Sp2] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                                [tempSd2,temp_psnr_Sd2] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
                                % Reference 3
                                sub_block_p = ref3(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                                sub_block_d = ref3(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
                                [tempSp3,temp_psnr_Sp3] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                                [tempSd3,temp_psnr_Sd3] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
                                % Reference 4
                                sub_block_p = ref4(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
                                sub_block_d = ref4(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
                                [tempSp4,temp_psnr_Sp4] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
                                [tempSd4,temp_psnr_Sd4] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
                                
                                tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
                                tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
                                temp_psnr_Sp = 1/4*(temp_psnr_Sp1 + temp_psnr_Sp2 + temp_psnr_Sp3 + temp_psnr_Sp4);
                                temp_psnr_Sd = 1/4*(temp_psnr_Sd1 + temp_psnr_Sd2 + temp_psnr_Sd3 + temp_psnr_Sd4);
                                
                                Sp_ref(ii,jj,:) = (tempSp(rang));
                                Sd_ref(ii,jj,:) = (tempSd(rang));
                                psnr_Sp_ref(ii,jj,:) = temp_psnr_Sp*ones(p,1);
                                psnr_Sd_ref(ii,jj,:) = temp_psnr_Sd*ones(p,1);

                                
%                                 env_data_block_sam  = env_rfdata_sam_roi(pixi,pixj,:);
%                                 SNR = mean(env_data_block_sam)/std(env_data_block_sam);
%                                 og_SNR(i,j) = SNR;
%                                 delta_SNR(i,j) = abs(SNR - 1.91).*100./1.91;
                                
                                if jj==floor(2.5*n/5)
                                    if ii==ceil(m/10)
                                        figure(103)
                                        %subplot(311);
                                        %set(gca,'FontSize',sl);
                                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'k--');   %
                                        title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                        axis([0 25 -70 0]);
                                        set(gca,'FontSize',font);
                                    elseif ii==round(m/2)
                                        figure(103);
                                        hold on;
                                        %subplot(312);
                                        set(gca,'FontSize',font);
                                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'r--');   %
                                        title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                        axis([0 25 -70 0]);
                                        set(gca,'FontSize',font);
                                    elseif ii==floor(9*m/10)
                                        figure(103);
                                        hold on;
                                        %subplot(313);
                                        set(gca,'FontSize',font);
                                        plot(band((1:NFFT/2+1))*1e-6,10*log10(tempSp((1:NFFT/2+1))/max(tempSp((1:NFFT/2+1)))),'b--');   %
                                        title('Spectrum REFERENCE'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
                                        axis([0 25 -70 0]);
                                        set(gca,'FontSize',font);
                                        %pause
                                        legend(' REF Top','REF Half','REF Bottom');
                                        hold off
                                    end
                                    
                                end
                                
                                
                            end
                        end
                        
                        diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;

                 figure(103); 
                 legend('Top','Half','Bottom'); title('SPECTRA')  
                 set(gca,'FontSize',20);

                % Au = b
                b = (log(Sp) - log(Sd)) - (diffraction_compensation);
                
                
                A1 = kron( 4*L*f , speye(m*n) );
                A2 = kron( ones(size(f)) , speye(m*n) );
                A = [A1 A2];
                

                SLD_term = b;

           SLD.SLD_term = SLD_term;
           SLD.band = f;
           SLD.A = A;
           SLD.A1 = A1;
           SLD.A2 = A2;
           SLD.x_ori = DATA.x;
           SLD.z_ori = DATA.z;
           SLD.z = z;
           SLD.x = x;
           SLD.x_ACS = x_ACS;
           SLD.z_ACS = z_ACS;
           SLD.zp_zd = -L;
           SLD.DATA_ROI = sam1;
           SLD.DATA_Bmode_ROI = my_RF2Bmode(sam1);
           
           SLD.DATA_Bmode = Im_db_ori;


            reshapedSLD = reshape(SLD.SLD_term, [], size(SLD.SLD_term, 3));
            meanValues = mean(reshapedSLD, 1, 'omitnan');
            stdValues = std(reshapedSLD, 0, 1, 'omitnan');
            
            ratios = abs(meanValues) ./ stdValues;
            SLD.SNR_ratios = ratios;
        
            figure, 
            plot(SLD.band', SLD.SNR_ratios, 'b.-'), title('Weights Julien'), grid, 
            xlabel('Frecuency [MHz]'), ylabel('SNR');
end


