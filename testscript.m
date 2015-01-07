
% PAPR curves for clipped SLM for varying S and number of candidates U
% and also for serial SLM + clipping 

%% --------------------------------------------------------------------------------------------------
% setup
% profile on


% close all
clear all
format short g
rand('state', sum(100*clock))  %#ok


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% waveform parameters + setup 
N_DFT = 256;            % FFT size
R     = 64;             % number of reserved tones
n_data_carriers             = N_DFT - R;
pow_signal                  = n_data_carriers/N_DFT;            % power of the OFDM signal
constellation_name          = '16QAM';
[ n_bits_per_carrier, M ]   = get_constellation_parameters(constellation_name);
Eb                          = 1/n_bits_per_carrier;
gamma_val                   = 1.4;      % value of clipping threshold gamma


% set up reserved tone pattern within OFDM symbol
rand('state', N_DFT) %#ok
tmp = randperm(N_DFT);
reserved_tone_positions = sort(tmp(1:R));
reserved_tone_pattern = false(1,N_DFT);
reserved_tone_pattern(reserved_tone_positions) = true;
rand('state', sum(100*clock)) %#ok



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% simulation parameters + setup
n_OFDM_symbols      = 1e2;%3e6;%2e6;    % 1e6 is approx. 3.5 hours computing time
S                   = 10;        % sparsity = number of clipped values
S_range             = 2:2:10;   % sparsity as a variable parameter
n_S                 = length(S_range);
n_bits              = n_data_carriers * n_bits_per_carrier;
n_SLM_runs          = 32;      % number of iterations of the dSLM procedure
% pre-calc rotation matrix for SLM
    rand('state',0) %#ok
    rot_mat = [   ones(N_DFT,1), ...      % 1st candidate is the original OFDM_symbol
                    1i.^randi(4, N_DFT, n_SLM_runs-1) ]; % axis-aligned QPSK
    rand('state', sum(100*clock))  %#ok

EbNo_dB_range      	= 0:2:16;
n_EbNo              = length(EbNo_dB_range);




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% generate signal

% init memory
PAPR_unclipped        = zeros(1,   n_OFDM_symbols);
idx_peaks             = zeros(S,   n_OFDM_symbols, 'uint16');
PAPR_ZERO             = zeros(1,   n_OFDM_symbols);
PAPR_MIN              = zeros(1,   n_OFDM_symbols);
PAPR_MAX              = zeros(1,   n_OFDM_symbols);
PAPR_SLM_2            = zeros(1,   n_OFDM_symbols);
PAPR_SLM_4            = zeros(1,   n_OFDM_symbols);
PAPR_SLM_8            = zeros(1,   n_OFDM_symbols);
PAPR_SLM_16           = zeros(1,   n_OFDM_symbols);
PAPR_SLM_32           = zeros(1,   n_OFDM_symbols);
PAPR_clSLM_2          = zeros(n_S, n_OFDM_symbols);
PAPR_clSLM_4          = zeros(n_S, n_OFDM_symbols);
PAPR_clSLM_8          = zeros(n_S, n_OFDM_symbols);
PAPR_clSLM_16         = zeros(n_S, n_OFDM_symbols);
PAPR_clSLM_32         = zeros(n_S, n_OFDM_symbols);
PAPR_serial_clSLM_2   = zeros(n_S, n_OFDM_symbols);
PAPR_serial_clSLM_4   = zeros(n_S, n_OFDM_symbols);
PAPR_serial_clSLM_8   = zeros(n_S, n_OFDM_symbols);
PAPR_serial_clSLM_16  = zeros(n_S, n_OFDM_symbols);
PAPR_serial_clSLM_32  = zeros(n_S, n_OFDM_symbols);


tic

for k=1:n_OFDM_symbols

    % create random bit sequence
    bitstream = ( rand(1, n_bits) < 0.5 );

    % create stream of data symbols
    mapped_signal = mapper(bitstream, constellation_name);


    % allocate to OFDM symbols
    %   columns = OFDM symbols; horiz. direction = time
    OFDM_signal_frq_domain = zeros(N_DFT, 1);
    OFDM_signal_frq_domain(~reserved_tone_pattern(:)) = mapped_signal(:);


    % time domain signal
    OFDM_signal_unclipped = ifft(OFDM_signal_frq_domain)*sqrt(N_DFT);

    
    % calc PAPR of unclipped OFDM signal
    PAPR_unclipped(k) = calc_PAPR(OFDM_signal_unclipped);    


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % clipping only
    [ ~, idx_OFDM_signal_sorted ] = sort(OFDM_signal_unclipped, 'descend');
    idx_peaks(1:S,k) = idx_OFDM_signal_sorted(1:S);     % remember positions of largest peaks
    
    
    %   Strategy ZERO - clipping to zeros
    OFDM_signal_clipped_ZERO = OFDM_signal_unclipped;
    OFDM_signal_clipped_ZERO(idx_peaks(:,k)) = 0;          % clip
    PAPR_ZERO(k) = calc_PAPR(OFDM_signal_clipped_ZERO);    % calc PAPR

    

    for idx_S=1:n_S
        
        S = S_range(idx_S);
        
            
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % SLM with clipping
        %   2 candidades
        [ ~, ~, PAPR_clSLM_2(idx_S,k), ~, ~ ] = clipped_SLM(OFDM_signal_frq_domain, rot_mat(:,1:2), S);
        %   4 candidades
        [ ~, ~, PAPR_clSLM_4(idx_S,k), ~, ~ ] = clipped_SLM(OFDM_signal_frq_domain, rot_mat(:,1:4), S);
        %   8 candidades
        [ ~, ~, PAPR_clSLM_8(idx_S,k), ~, ~ ] = clipped_SLM(OFDM_signal_frq_domain, rot_mat(:,1:8), S);
        %   16 candidades
        [ ~, ~, PAPR_clSLM_16(idx_S,k), ~, ~ ] = clipped_SLM(OFDM_signal_frq_domain, rot_mat(:,1:16), S);
        %   32 candidades
        [ ~, ~, PAPR_clSLM_32(idx_S,k), ~, ~ ] = clipped_SLM(OFDM_signal_frq_domain, rot_mat(:,1:32), S);


    
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % serial SLM plus clipping
        %   2 candidades
        [ output_OFDM_symbol, ~, ~ ]  = basic_SLM(OFDM_signal_frq_domain, rot_mat(:,1:2));
        output_OFDM_symbol = ifft(output_OFDM_symbol) * sqrt(N_DFT);
        tmp = sort(real(output_OFDM_symbol).^2 + imag(output_OFDM_symbol).^2, 'descend');
        PAPR_serial_clSLM_2(idx_S,k) = tmp(S+1) / sum(tmp(S+1:end)) * N_DFT;
        
        %   4 candidades
        [ output_OFDM_symbol, ~, ~ ]  = basic_SLM(OFDM_signal_frq_domain, rot_mat(:,1:4));
        output_OFDM_symbol = ifft(output_OFDM_symbol) * sqrt(N_DFT);
        tmp = sort(real(output_OFDM_symbol).^2 + imag(output_OFDM_symbol).^2, 'descend');
        PAPR_serial_clSLM_4(idx_S,k) = tmp(S+1) / sum(tmp(S+1:end)) * N_DFT;

        %   8 candidades
        [ output_OFDM_symbol, ~, ~ ]  = basic_SLM(OFDM_signal_frq_domain, rot_mat(:,1:8));
        output_OFDM_symbol = ifft(output_OFDM_symbol) * sqrt(N_DFT);
        tmp = sort(real(output_OFDM_symbol).^2 + imag(output_OFDM_symbol).^2, 'descend');
        PAPR_serial_clSLM_8(idx_S,k) = tmp(S+1) / sum(tmp(S+1:end)) * N_DFT;

        %   16 candidades
        [ output_OFDM_symbol, ~, ~ ]  = basic_SLM(OFDM_signal_frq_domain, rot_mat(:,1:16));
        output_OFDM_symbol = ifft(output_OFDM_symbol) * sqrt(N_DFT);
        tmp = sort(real(output_OFDM_symbol).^2 + imag(output_OFDM_symbol).^2, 'descend');
        PAPR_serial_clSLM_16(idx_S,k) = tmp(S+1) / sum(tmp(S+1:end)) * N_DFT;

        %   32 candidades
        [ output_OFDM_symbol, ~, ~ ]  = basic_SLM(OFDM_signal_frq_domain, rot_mat(:,1:32));
        output_OFDM_symbol = ifft(output_OFDM_symbol) * sqrt(N_DFT);
        tmp = sort(real(output_OFDM_symbol).^2 + imag(output_OFDM_symbol).^2, 'descend');
        PAPR_serial_clSLM_32(idx_S,k) = tmp(S+1) / sum(tmp(S+1:end)) * N_DFT;
    
    
    end
    
    
    
    
end

toc








save results_script_07




%%


load results_script_07


% profile on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  plot histograms/ccdf
axis_PAPR = 4:.1:14;


%  observed PAPRs
ccdf_PAPR_clSLM_2  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_clSLM_4  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_clSLM_8  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_clSLM_16 = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_clSLM_32 = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_serial_clSLM_2  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_serial_clSLM_4  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_serial_clSLM_8  = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_serial_clSLM_16 = zeros(n_S,numel(axis_PAPR));
ccdf_PAPR_serial_clSLM_32 = zeros(n_S,numel(axis_PAPR));


for idx_S=1:n_S
    
    
    ccdf_PAPR_clSLM_2(idx_S,:)      = histc(10*log10(PAPR_clSLM_2(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_clSLM_2(idx_S,:)    = ccdf_PAPR_clSLM_2(idx_S,end:-1:1);
      ccdf_PAPR_clSLM_2(idx_S,:)    = cumsum(ccdf_PAPR_clSLM_2(idx_S,:) );
      ccdf_PAPR_clSLM_2(idx_S,:)    = ccdf_PAPR_clSLM_2(idx_S,end:-1:1);

    ccdf_PAPR_clSLM_4(idx_S,:)      = histc(10*log10(PAPR_clSLM_4(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_clSLM_4(idx_S,:)    = ccdf_PAPR_clSLM_4(idx_S,end:-1:1);
      ccdf_PAPR_clSLM_4(idx_S,:)    = cumsum(ccdf_PAPR_clSLM_4(idx_S,:) );
      ccdf_PAPR_clSLM_4(idx_S,:)    = ccdf_PAPR_clSLM_4(idx_S,end:-1:1);

    ccdf_PAPR_clSLM_8(idx_S,:)      = histc(10*log10(PAPR_clSLM_8(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_clSLM_8(idx_S,:)    = ccdf_PAPR_clSLM_8(idx_S,end:-1:1);
      ccdf_PAPR_clSLM_8(idx_S,:)    = cumsum(ccdf_PAPR_clSLM_8(idx_S,:) );
      ccdf_PAPR_clSLM_8(idx_S,:)    = ccdf_PAPR_clSLM_8(idx_S,end:-1:1);

    ccdf_PAPR_clSLM_16(idx_S,:)     = histc(10*log10(PAPR_clSLM_16(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_clSLM_16(idx_S,:)   = ccdf_PAPR_clSLM_16(idx_S,end:-1:1);
      ccdf_PAPR_clSLM_16(idx_S,:)   = cumsum(ccdf_PAPR_clSLM_16(idx_S,:) );
      ccdf_PAPR_clSLM_16(idx_S,:)   = ccdf_PAPR_clSLM_16(idx_S,end:-1:1);

    ccdf_PAPR_clSLM_32(idx_S,:)     = histc(10*log10(PAPR_clSLM_32(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_clSLM_32(idx_S,:)   = ccdf_PAPR_clSLM_32(idx_S,end:-1:1);
      ccdf_PAPR_clSLM_32(idx_S,:)   = cumsum(ccdf_PAPR_clSLM_32(idx_S,:) );
      ccdf_PAPR_clSLM_32(idx_S,:)   = ccdf_PAPR_clSLM_32(idx_S,end:-1:1);

    
    ccdf_PAPR_serial_clSLM_2(idx_S,:)      = histc(10*log10(PAPR_serial_clSLM_2(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_serial_clSLM_2(idx_S,:)    = ccdf_PAPR_serial_clSLM_2(idx_S,end:-1:1);
      ccdf_PAPR_serial_clSLM_2(idx_S,:)    = cumsum(ccdf_PAPR_serial_clSLM_2(idx_S,:) );
      ccdf_PAPR_serial_clSLM_2(idx_S,:)    = ccdf_PAPR_serial_clSLM_2(idx_S,end:-1:1);

    ccdf_PAPR_serial_clSLM_4(idx_S,:)      = histc(10*log10(PAPR_serial_clSLM_4(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_serial_clSLM_4(idx_S,:)    = ccdf_PAPR_serial_clSLM_4(idx_S,end:-1:1);
      ccdf_PAPR_serial_clSLM_4(idx_S,:)    = cumsum(ccdf_PAPR_serial_clSLM_4(idx_S,:) );
      ccdf_PAPR_serial_clSLM_4(idx_S,:)    = ccdf_PAPR_serial_clSLM_4(idx_S,end:-1:1);

    ccdf_PAPR_serial_clSLM_8(idx_S,:)      = histc(10*log10(PAPR_serial_clSLM_8(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_serial_clSLM_8(idx_S,:)    = ccdf_PAPR_serial_clSLM_8(idx_S,end:-1:1);
      ccdf_PAPR_serial_clSLM_8(idx_S,:)    = cumsum(ccdf_PAPR_serial_clSLM_8(idx_S,:) );
      ccdf_PAPR_serial_clSLM_8(idx_S,:)    = ccdf_PAPR_serial_clSLM_8(idx_S,end:-1:1);

    ccdf_PAPR_serial_clSLM_16(idx_S,:)     = histc(10*log10(PAPR_serial_clSLM_16(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_serial_clSLM_16(idx_S,:)   = ccdf_PAPR_serial_clSLM_16(idx_S,end:-1:1);
      ccdf_PAPR_serial_clSLM_16(idx_S,:)   = cumsum(ccdf_PAPR_serial_clSLM_16(idx_S,:) );
      ccdf_PAPR_serial_clSLM_16(idx_S,:)   = ccdf_PAPR_serial_clSLM_16(idx_S,end:-1:1);

    ccdf_PAPR_serial_clSLM_32(idx_S,:)     = histc(10*log10(PAPR_serial_clSLM_32(idx_S,:)),       axis_PAPR)/n_OFDM_symbols;
      ccdf_PAPR_serial_clSLM_32(idx_S,:)   = ccdf_PAPR_serial_clSLM_32(idx_S,end:-1:1);
      ccdf_PAPR_serial_clSLM_32(idx_S,:)   = cumsum(ccdf_PAPR_serial_clSLM_32(idx_S,:) );
      ccdf_PAPR_serial_clSLM_32(idx_S,:)   = ccdf_PAPR_serial_clSLM_32(idx_S,end:-1:1);

end


% calc theoretical PAPR of Gaussian distributed frames
p0 = 10.^(axis_PAPR/10);
ccdf_PAPR_theoretical = 1-(1-exp(-p0)).^N_DFT;




%%


darkgreen = 0.6 * [ 0 1 0 ];
darkcyan  = 0.7 * [ 0 1 1 ];

% plot PAPR for clSLM
figure(41)
clf
p{1} = semilogy(   axis_PAPR, ccdf_PAPR_theoretical, '--k', 'LineWidth', 1.5) ;
hold on
p{2} =  semilogy(   axis_PAPR, ccdf_PAPR_clSLM_2(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_clSLM_4(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_clSLM_8(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_clSLM_16(3,:), 'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_clSLM_32(3,:), 'b', 'LineWidth', 1.5);
p{3} =  semilogy(   axis_PAPR, ccdf_PAPR_clSLM_16,  '--', 'LineWidth', 1.5, 'Color', darkgreen);
grid on
zoom on
axis([4,10,1e-4,1])
% axis([4,8,1e-4,1])
set(gcf, 'Position', [ 10 650 600 350 ])
xlabel('$10 \log_{10} ($PAPR$_0)$', 'interpreter', 'latex')
ylabel('ccdf$($PAPR$_0)$', 'interpreter', 'latex')


legend( [p{1}(1),p{2}(1),p{3}(1)] ,...
           'OFDM  unclipped    {}', ...
        [  '                      ';
           'clSLM(6,{\itU})       ';
           '  {\itU} = 2,4,8,16,32';
           '   (right to left)    '], ...
        [  '                      ';
           'clSLM({\itS},16)      ';
           '  {\itS} = 2,4,6,8,10 ';
           '   (right to left)    '], ...
        'location', 'southeast')
shg



% plot PAPR for serial SLM + clipping
figure(42)
clf
p{1} = semilogy(   axis_PAPR, ccdf_PAPR_theoretical, '--k', 'LineWidth', 1.5) ;
hold on
p{2} =  semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_2(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_4(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_8(3,:),  'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_16(3,:), 'b', 'LineWidth', 1.5);
        semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_32(3,:), 'b', 'LineWidth', 1.5);
p{3} =  semilogy(   axis_PAPR, ccdf_PAPR_serial_clSLM_16,  '--', 'LineWidth', 1.5, 'Color', darkgreen);
grid on
zoom on

axis([4,10,1e-4,1])
% axis([4,8,1e-4,1])
set(gcf, 'Position', [ 10 200 600 350 ])
xlabel('$10 \log_{10} ($PAPR$_0)$', 'interpreter', 'latex')
ylabel('ccdf$($PAPR$_0)$', 'interpreter', 'latex')


h_legend=legend( [p{1}(1),p{2}(1),p{3}(1)] ,...
           'OFDM  unclipped          {}', ...
        [  '                        ';
           'SLM + clipping (serial) ';
           '  {\itS} = 6            ';
           '  {\itU} = 2,4,8,16,32  ';
           '   (right to left)      '], ...
        [  '                        ';
           'SLM + clipping (serial) ';
           '  {\itU} = 16           ';
           '  {\itS} = 2,4,6,8,10   ';
           '   (right to left)      '], ...
        'location', 'southeast');

set(h_legend,'FontSize',9);

shg












