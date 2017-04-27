% Federal University of Rio de Janeiro - UFRJ
% Electrical Engineering Program - COPPE
% Signals, Multimedia, and Telecommunications Laboratory - SMT
%
% Author: Rafael da Silva Chaves
% email: rafael.chves@smt.ufrj.br
%
% Advisors: Wallace Alves Martins and Paulo S. R. Diniz
% email: {wallace.martins, diniz}@smt.ufrj.br
%
% Abstract: This function generates a impulse response for underwater
% acoustic (UWA) channel, taking into account Doppler effect. This function
% uses the UWA channel model presented in [1], and use the metodolgy
% presented in [2], which is based in the UWA channels models described in
% [1].
%
% References:
%
% [1] - S. Zhou and Z. Wang, OFDM for Underwater Acoustic Communications.
% Chichester: John Wiley & Sons, May 2014.
%
% [2] - R. S. Chaves, "Modeling and Simulation of Underwater Acoustic
% Communication Systems (in Portuguese)", Graduation Thesis, Federal
% Universty of Rio de Janeiro, Rio de Janeiro, 2016.

function [h,delay_bar,gain_tap,varargout] = acousticChannel(setup,gain,delay,doppler)
% [h,delay_bar,gain_tap,varargout] = acousticChannel(setup,gain,delay,doppler)
% 
% Inputs:
%        --'setup' is a structure with the parameters for channel.
%
%        'setup.Ts' is the field with sampling rate of UWA channel, must be
%        a real value.
%        'setup.paths' is the field with the number of paths of UWA
%        channel, must be a positive integer value.
%        'setup.delayspread' is the field with the delay spread of UWA
%        channel, muste be a real value.
%
%        --'gain' is a structure with the parameters for gains.
%
%        'gain.attenuation' is the field with the attenuation that occur
%        during delay spread, must be a real value in dB.
%
%        --'delay' is a structure with parameters for delays.
%
%        'delay.mean' is the field with the mean of the exponential pdf for
%        delays, must be a real value.
%
%        --'doppler' is a structure with parameters for doppler.
%
%        'doppler.type' is the field for the Doppler effect in channel,
%        must be a string with the following values:
%              'none' - for an UWA channel with none Doppler scaling 
%                       factor (DSF).
%              'uniform' - for an UWA channel with uniform DSF.
%              'non-uniform' - for an UWA channel with non-uniform DSF.
%        'doppler.velocity' is the field for the relative motion between
%        transmitter and receptor, must be a real value for uniform DSF,
%        and must be an array with the same size of multipaths for 
%        non-uniform DSF.
%
% Outputs:
%         --'h' is a real array with impulse response of UWA channel.
%
%         --'delay_bar' is a real array with the delay for each multipath.
%
%         --'gain_tap' is a real array with the gains for each multipath.
%
%         --'varargout':
%
%         varargout{1} = Q is the downsamplig factor of the DSF.
%         varargout{2} = M is the upsampling factor of the DSF.

% Generating the distrubition of delays

distribution_delay = makedist('exponential','mu',delay.mean);  % Generating an exponential pdf for Dtau

Dtau        = random(distribution_delay,setup.paths,1);        % Vector of Dtau exponentially distributed 
Dtau_index  = ceil(Dtau/setup.Ts);                             % Indexes of Dtau vector in disrete time

delay_index = cumsum(Dtau_index);                              % Calculating indexes of multipaths for discrete time
delay       = (delay_index - 1)*setup.Ts;                      % Calculating delays of multipaths

% Generating the distribution of gains

alpha = log(10^(gain.attenuation/10))/setup.delayspread;       % Exponential attenuation power factor

gain_variance = exp(-alpha*delay);                             % Calculating the gain ariance for each multipath                      
gain_tap      = raylrnd(sqrt(gain_variance*2/(4 - pi)));       % Calculating the gain for each multipath

gain_tap = gain_tap/norm(gain_tap);                            % Normalizing gain taps

% Doppler effect

c = 1500; % Sound speed

% Testing type of doppler effect in channel

if(strcmp(doppler.type,'none'))
    a_max = 0;
    varargout{1} = 1;
    varargout{2} = 1;
elseif(strcmp(doppler.type,'uniform'))
    a_max = doppler.velocity/c;
    [Q,M] = rat(1+a_max);
    varargout{1} = Q;
    varargout{2} = M;
elseif (strcmp(type,'non-uniform'))
    a = doppler.velocity/c;
    a_max = max(a);
    [Q,M] = rat(1+a_max);
    varargout{1} = Q;
    varargout{2} = M;
else
    error('Invalid type.')
end

delay_index_bar = ceil(delay_index./(1 + a_max));              % Calculating new delays indexes corrected by DSF
delay_bar       = (delay_index_bar - 1)*setup.Ts;              % Calculating new delays corrected by DSF

h(delay_index_bar) = gain_tap;                                 % Channel impulse response

end