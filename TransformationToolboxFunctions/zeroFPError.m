function H = zeroFPError(H,ZERO)
%zeroFPError zeros round-off terms caused by floating point errors
%   zeroFPError(H) simplifies a matrix H by zeroing floating point errors 
%   with magnitudes less than the 10*eps(class(H)) for single and double 
%   precision variables, and a default value of 1e-10 otherwise.
%
%   zeroFPError(H,ZERO) simplifies a matrix H by zeroing floating point  
%   errors with magnitudes less than the value specified in ZERO.
%
%   See also sym/simplify sym/vpa eps isZero
%
%   M. Kutzer 27Oct2014, USNA

%Updates
%   11Sep2014 - Remove extra VPA call.
%   22Sep2014 - Updated to check for use of symbolic variables prior to 
%               calling SIMPLIFY.
%   22Sep2014 - Added SIMPLIFY call prior to calling VPA to ensure trig 
%               simplification.
%   09Oct2014 - Updated to accommodate inputs that are fractions.
%   27Oct2014 - Round coefficients based on zero. This includes additional 
%               VPA for decimal coefficients in output.
%   27Oct2014 - Added "See also".
%   13May2015 - Updated definition of zero and "See also".
%   07Nov2018 - Updated documentation of default zero.

%% Set Defaults
if nargin < 2
    if strcmpi(class(H),'double') || strcmpi(class(H),'single')
        ZERO = 10*eps(class(H));
    else
        ZERO = 1e-10; % assume magnitudes smaller than this are zero
    end
end 

%% Zero floating point error
if strcmpi( class(H), 'sym') % apply for symbolic inputs
    % Simplify initial input (prior to VPA call)
    H = simplify(H);

    % Numerically evaluate terms
    dgts = ceil( -log10(ZERO) ); % specify digits based on ZERO
    H = vpa(H,dgts+1);           % numerically evaluate input
    
    % Zero
    for i = 1:numel(H)
        h = H(i);
        [n,d] = numden(h); % account for expressions that are fractions
        
        % Numerator
        [c,s] = coeffs(n);
        if ~isempty(c) && ~isempty(s)
            for k = 1:numel(c)
                if abs( c(k) ) < ZERO
                    c(k) = 0;
                else % round based on zero
                    c(k) = round(c(k)*(1*10^dgts))*(1*10^-dgts);
                end
            end
            n = s*c';
        else
            n = 0;
        end
        
        % Denominator 
        if d ~= 1
            [c,s] = coeffs(d);
            for k = 1:numel(c)
                if abs( c(k) ) < ZERO
                    c(k) = 0;
                else % round based on zero
                    c(k) = round(c(k)*(1*10^dgts))*(1*10^-dgts);
                end
            end
            d = s*c';
            
            % Combine
            H(i) = simplify(n/d);
        else % for denominator of 1
            H(i) = n;
        end

    end
    H = vpa(H);
else % apply for double and single precision
    for i = 1:numel(H)
        h = H(i);
        % Zero
        if abs(h) < ZERO
            H(i) = 0;
        end
    end
end