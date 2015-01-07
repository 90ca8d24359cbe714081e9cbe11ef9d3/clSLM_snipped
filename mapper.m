function mapped_signal = mapper(input_stream, constellation_name)
% BPSK/4-QAM/16-QAM mapping


const_name = lower(constellation_name);

switch const_name
    case {'bpsk'}
        
        mapped_signal = 1 - 2*input_stream ;
                
    case {'4qam', '4-qam', 'qpsk'}
        
        mapped_signal = (1 - 2*input_stream(1:2:end) + 1i * (1 - 2*input_stream(2:2:end))) /sqrt(2);
        
    case {'16qam', '16-qam'}
        
        mapped_signal = ( ...
                   (1 - 2*input_stream(1:4:end)) .* (3 - 2*input_stream(3:4:end)) +     ...
                   (1 - 2*input_stream(2:4:end)) .* (3 - 2*input_stream(4:4:end)) * 1i  ...
                        ) / sqrt(10);  

    otherwise 
        error('Invalid mapping type specified')
       
end


end


