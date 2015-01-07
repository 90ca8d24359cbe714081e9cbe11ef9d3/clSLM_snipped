function [ output_OFDM_symbol, idx_candidate, PAPR_min ] = ...
                            basic_SLM(OFDM_symbol, n_SLM_runs__or__rot_mat)
% Basic version of SLM (selected mapping)
% 
% 
% usage: basic_SLM(OFDM_symbol, n_SLM_runs__or__rot_mat)
%    input: OFDM_symbol             - a column vector(!!) of length N_DFT in freq. domain 
%           n_SLM_runs__or__rot_mat - either a scalar (number of candidates)
%                                     or a matrix of size  N_DFT x number of candidates
% 
% 
%   output: output_OFDM_symbol      - selected candidate in freq. domain 
%           idx_candidate           - index/number of the selected candidate
%           PAPR_min                - PAPR of the selected candidate reached with SLM
% 
% 



% init stuff
N_DFT = size(OFDM_symbol,1); % OFDM symbol has to be a column vector!

if isscalar(n_SLM_runs__or__rot_mat)
    n_SLM_runs = n_SLM_runs__or__rot_mat;        % input = number of candidates
    % create reproduceable rotation matrix 
    rand('state',0) %#ok
    rotations = [   ones(N_DFT,1), ...      % 1st candidate is the original OFDM_symbol
                    1i.^randi(4, N_DFT, n_SLM_runs-1) ]; % axis-aligned QPSK
else
    rotations = n_SLM_runs__or__rot_mat;         % input = matrix of rotations
    n_SLM_runs = size(rotations,2);
end


% rotate OFDM symbol
% rotated_OFDM_symbols = rotations .* OFDM_symbol(:, ones(1,n_SLM_runs));
rotated_OFDM_symbols = bsxfun(@times, rotations, OFDM_symbol);


% % calc PAPRs in time domain
% time_domain_OFDM_symbols = ifft(rotated_OFDM_symbols);
% 
% PAPR_results = zeros(1,n_SLM_runs);
% for k=1:n_SLM_runs
%     PAPR_results(k) = calc_PAPR(time_domain_OFDM_symbols(:,k));
% end


% calc PAPRs in time domain
tmp = ifft(rotated_OFDM_symbols);
tmp = real(tmp).^2 + imag(tmp).^2;  % |data|²

PAPR_results(1,:) = max(tmp) ./ sum(tmp) * N_DFT; % PAPR all





% find minimal PAPR and output respective OFDM symbol
[ PAPR_min, idx_candidate ]  = min(PAPR_results);
output_OFDM_symbol           = rotated_OFDM_symbols(:,idx_candidate);








end