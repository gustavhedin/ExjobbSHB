classdef ADRev <  handle
    % This is a class defineing ADRev objects, which are used to calculate
    % function value and derivatives using the reverso mode of Automatic
    % differentiation. 
    % 
    % Implemented by Gustav Hedin, 2018. With inspiration from Laksh Gupta.
    % 
    % Methods in this class are: 
    %   addition          adr_sub()
    %   subtraction       adr_sub()
    %   multiplication    adr_mul()
    %   division          adr_div()
    %   natural exponent  adr_exp() 
    %   natural logarithm adr_ln()  
    %   ...
    %

    properties
        value
        derivative
        derivativeOp
        parents
    end
    
    methods
        % Constructor: 
        function obj = ADRev(val,der)
            if nargin == 1
                obj.value = val;
                obj.derivative = 0;
                obj.derivativeOp = @ adr_constD;
                obj.parents = [];
            else
                obj.value = val;
                obj.derivative = der;
                obj.derivativeOp = @ adr_constD;
                obj.parents = [];
            end
        end
        
        function modified_parents = adr_constD(~, ~)
            %adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            %adNodes(2).derivative = adNodes(2).derivative + prevDerivative;
            modified_parents = [];
        end
        
        % Fuctions for defineing the basic calculation operations for
        % ADRev-objects:
        
        % Addition
        function result = adr_add(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value + y.value);
            result.derivativeOp = @ adr_addD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_addD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative;
            modified_parents = adNodes;
        end
        
        % Subtraction
        function result = adr_sub(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value - y.value);
            result.derivativeOp = @ adr_subD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_subD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative - prevDerivative;
            modified_parents = adNodes;
        end
        
        % Multiplication:
        function result = adr_mul(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value .* y.value);
            result.derivativeOp = @ adr_mulD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_mulD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*adNodes(2).value;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative.*adNodes(1).value;
            modified_parents = adNodes;
        end
        
        % Division
        function result = adr_div(x,y)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            elseif ~isa(y,'ADRev')
                y = ADRev(y);
            end
            result = ADRev(x.value ./ y.value);
            result.derivativeOp = @ adr_divD;
            result.parents = [x y];
        end
        
        function modified_parents = adr_divD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative*(1/adNodes(2).value);
            adNodes(2).derivative = adNodes(2).derivative + ...
                prevDerivative.*(-adNodes(1).value./adNodes(2).value.^2);
            modified_parents = adNodes;
        end
        
        % Natural exponent
        function result = adr_exp(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(exp(x.value));
            result.derivativeOp = @ adr_expD;
            result.parents = x;
        end
        
        function modified_parents = adr_expD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*exp(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % Natural logarithm
        function result = adr_ln(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(log(x.value));
            result.derivativeOp = @ adr_lnD;
            result.parents = x;
        end
        
        function modified_parents = adr_lnD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*(1/adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % Power
        function result = adr_pow(x,p)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(x.value.^p);
            switch p
                case 2
                    result.derivativeOp = @ adr_powD_2;
                case 3
                    result.derivativeOp = @ adr_powD_3;
                case 4
                    result.derivativeOp = @ adr_powD_4;
                otherwise
                    error('x^p for p<4 not implemented yet. use a lower p or extend implementation in ADRev')
            end
            result.parents = x;
        end
        
        function modified_parents = adr_powD_2(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 2*prevDerivative.*adNodes(1).value;
            modified_parents = adNodes;
        end
        
        function modified_parents = adr_powD_3(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 3*prevDerivative.*adNodes(1).value^2;
            modified_parents = adNodes;
        end
        
        function modified_parents = adr_powD_4(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + 4*prevDerivative.*adNodes(1).value^3;
            modified_parents = adNodes;
        end
        
        % sine
        function result = adr_sin(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(sin(x.value));
            result.derivativeOp = @ adr_sinD;
            result.parents = x;
        end
        
        function modified_parents = adr_sinD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*cos(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % cosine
        function result = adr_cos(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(cos(x.value));
            result.derivativeOp = @ adr_cosD;
            result.parents = x;
        end
        
        function modified_parents = adr_cosD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*sin(adNodes(1).value);
            modified_parents = adNodes;
        end
        
        % square root
        function result = adr_sqrt(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            result = ADRev(sqrt(x.value));
            result.derivativeOp = @ adr_sqrtD;
            result.parents = x;
        end
        
        function modified_parents = adr_sqrtD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*1/(2*sqrt(adNodes(1).value));
            modified_parents = adNodes;
        end
        
        % step-function
        function result = adr_step(x)
            if ~isa(x,'ADRev')
                x = ADRev(x);
            end
            v = x.value;
            unistep = v>=0;
            result = ADRev(unistep);
            result.derivativeOp = @ adr_stepD;
            result.parents = x;
        end
        
        function modified_parents = adr_stepD(prevDerivative, adNodes)
            a = 0.001; % S?tta denna parameter n?gon annanstans? i properties?
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative.*(exp(-adNodes(1).value.^2/(2*a))/(sqrt(a*pi)));
            modified_parents = adNodes;
        end
        
    end
   
end

