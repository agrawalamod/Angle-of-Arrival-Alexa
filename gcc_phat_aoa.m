
function [deg1, deg2] = gcc_phat_aoa(fy, Fs, c, res, azi_2d, ele_2d, mic_dis)
    fy_tau1 = [];
    fy_tau2 = [];
    for ref_mic = 1:6
        refsig = fy(:,ref_mic);
        for idx=1:6
            if(ref_mic == idx)
                tau1=0;
                tau2=0;
            else
                %tau1 = gccphat(refsig, fy(:,idx), Fs);
                tau2 = interpolated_gccphat(refsig, fy(:,idx), Fs);
                tau1=tau2;
            end
            fy_tau1 = [fy_tau1; tau1];
            fy_tau2 = [fy_tau2; tau2];
        end
    end
    fy_dis1 = fy_tau1*c;
    fy_dis2 = fy_tau2*c;
    
    scores1 = zeros(size(azi_2d));
    scores2 = zeros(size(azi_2d));

    
    for A=1:size(azi_2d,2)
        for E=1:size(azi_2d,1)
            azi = azi_2d(E,A)*pi/180;
            ele = ele_2d(E,A)*pi/180;
            vec = [cos(azi)*cos(ele) sin(azi)*cos(ele) sin(ele)];
            scores1(E,A)= 1/norm(mic_dis*vec' - fy_dis1);
            scores2(E,A)= 1/norm(mic_dis*vec' - fy_dis2);            
        end
    end
    
    result1 = sum(scores1,1);
    result2 = sum(scores2,1);
    
    [max_val1, max_indx1] = max(result1);
    [max_val2, max_indx2] = max(result2);
    
    deg1 = res*max_indx1;
    deg2 = res*max_indx2;
    
end

function [tau] = interpolated_gccphat(y1, y2, Fs)
    esp = 0.000001;
    N = length(y1);
    
    ffta_y1 = fft(y1,2*N);
    fft_y1 = ffta_y1(1:(floor(length(ffta_y1)/2)+1));
    ffta_y2 = fft(y2,2*N);
    fft_y2 = ffta_y2(1:(floor(length(ffta_y2)/2)+1));
    
    fft_prod = fft_y1 .* conj(fft_y2);
    fft_prod_mag = abs(fft_prod);
    Gphat = fft_prod./(fft_prod_mag+esp);
    Rphat = irfft(Gphat, 1, 2*N);
    
    x = horzcat(Rphat(N+1:end),Rphat(1:N));
    lag = -N:N;
    %figure;
    %plot(1:size(x,2),x);
    %title('IFFT of GCC PHAT: Correlation Values')
    
   %splining
   [corr_val, corr_index] = max(x);
   tau_approx = lag(corr_index)/Fs;
   
   dk = 0.5*(x(corr_index-1) - x(corr_index+1)) / (x(corr_index+1) + x(corr_index-1) - 2*x(corr_index));
   if(abs(dk)>=1)
       dk=0;
   end
   %{
   if(tau_approx==0)
       dk=0;
   end
   %}
   t_offset = dk/Fs;
   tau = tau_approx + t_offset;
end

function [first_mic, mic_delays] = first_mic_solver(y)
    mic_delays = [];
    for i=1:6
        [res, lags] = xcorr(y(:,1), y(:,i));
        [max_v, max_i] = max(res);
        mic_delays = [mic_delays lags(max_i)];
    end
    [least_lag] = max(mic_delays);
    first_mic = find(mic_delays == least_lag);
end

function [irfft_result] = irfft(x,even, N)
     n = 0; % the output length
     s = 0; % the variable that will hold the index of the highest
            % frequency below N/2, s = floor((n+1)/2)
     if(even==1)
        n = 2 * (length(x) - 1 );
        s = length(x) - 1;
     else
        n = 2 * (length(x) - 1 )+1;
        s = length(x);
     end
     xn = zeros(1,n);
     xn(1:length(x)) = x;
     xn(length(x)+1:n) = conj(x(s:-1:2));
     irfft_result  = ifft(xn, N);
end
