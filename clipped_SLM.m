function [ output_OFDM_symbol, idx_candidate, PAPR_clipped, idx_clip, PAPR_unclipped, sorted_powers ] = ...
                            clipped_SLM(OFDM_symbol, n_SLM_runs__or__rot_mat, S)
% Clipping version of SLM (selected mapping), allowing S peaks in time
% domain, which are not considered during PAPR calculation/optimization
% Note: only ZERO clipping is considered, MAX and MIN can be obtained from
%       the result
% 
% usage: basic_SLM(OFDM_symbol, n_SLM_runs__or__rot_mat, S)
%    input: OFDM_symbol             - a column vector(!!) of length N_DFT in freq. domain 
%           n_SLM_runs__or__rot_mat - either a scalar (number of candidates)
%                                     or a matrix of size  N_DFT x number of candidates
%           S                       - number of clipped peaks
% 
%   output: output_OFDM_symbol      - selected candidate in freq. domain (unclipped)
%           idx_candidate           - index/number of the selected candidate
%           PAPR_clipped            - PAPR of the selected candidate reached with clipped SLM
%           idx_clip                - indices of the clipped peaks in time domain
%           PAPR_unclipped          - PAPR of the selected candidate if not clipped 
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
rotated_OFDM_symbols = rotations .* OFDM_symbol(:, ones(1,n_SLM_runs));


% % calc PAPRs in time domain
% time_domain_OFDM_symbols = ifft(rotated_OFDM_symbols);
% [ ~, idx ] = sort(abs(time_domain_OFDM_symbols), 'descend');
% idx = idx(1:S,:);
% 
% PAPR_results = zeros(2,n_SLM_runs);
% for k=1:n_SLM_runs
%     tmp = time_domain_OFDM_symbols(:,k);
%     PAPR_results(1,k) = calc_PAPR(tmp);     % PAPR unclipped
%     tmp(idx(:,k)) = 0;                      % clip to zero
%     PAPR_results(2,k) = calc_PAPR(tmp);     % PAPR clipped
% end



% calc PAPRs in time domain
tmp = ifft(rotated_OFDM_symbols);
tmp = real(tmp).^2 + imag(tmp).^2;  % |data|²
[ sorted_powers, idx ] = sort(tmp, 'descend');
idx = idx(1:S,:);

sum_clipped  = sum(sorted_powers(S+1:end,:));
PAPR_results = zeros(2,n_SLM_runs);
PAPR_results(1,:) = sorted_powers(1,:)   ./ (sum_clipped + sum(sorted_powers(1:S,:))); % PAPR all
PAPR_results(2,:) = sorted_powers(S+1,:) ./  sum_clipped ;                   % PAPR clipped
PAPR_results = PAPR_results * N_DFT;



% find minimal PAPR and output respective OFDM symbol
[ PAPR_clipped, idx_candidate ]  = min(PAPR_results(2,:));
output_OFDM_symbol               = rotated_OFDM_symbols(:,idx_candidate);

PAPR_unclipped = PAPR_results(1,idx_candidate);
idx_clip = idx(:, idx_candidate);


sorted_powers = sorted_powers * N_DFT;


end