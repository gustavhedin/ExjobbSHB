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
    %   natural exponent  adr_exp() Not yet implemented
    %   natural logarithm adr_ln()  Not yet implemented
    %   ...
    %

    properties
        value
        derivative
        derivativeOp
        parents
    end
    methods
        function obj = ADRev(val,der)
            if nargin == 1
                obj.value = val;
                obj.derivative = 0;
                obj.derivativeOp = @ ad_constD;
                obj.parents = [];
            else
                obj.value = val;
                obj.derivative = der;
                obj.derivativeOp = @ ad_constD;
                obj.parents = [];
            end
        end
        
        function output = ad_constD(prevDerivatives, adNodes)
            output = 0;
        end
        
        % Addition
        function result = adr_add(x,y)
            result = ADRev(x.value + y.value);
            result.derivativeOp = @ adr_addD;
            result.parents = [result.parents x y];
        end
        
        function modified_parents = adr_addD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative;
            modified_parents = adNodes;
        end
        
        % Subtraction
        function result = adr_sub(x,y)
            result = ADRev(x.value - y.value);
            result.derivativeOp = @ adr_subD;
            result.parents = [result.parents x y];
        end
        
        function modified_parents = adr_subD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative;
            adNodes(2).derivative = adNodes(2).derivative - prevDerivative;
            modified_parents = adNodes;
        end
        
        % Multiplication:
        function result = adr_mul(x,y)
            result = ADRev(x.value * y.value);
            result.derivativeOp = @ adr_mulD;
            result.parents = [result.parents x y]; % R?tt?
        end
        
        function modified_parents = adr_mulD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative*adNodes(2).value;
            adNodes(2).derivative = adNodes(2).derivative + prevDerivative*adNodes(1).value;
            modified_parents = adNodes;
        end
        
        % Division
        function result = adr_div(x,y)
            result = ADRev(x.value / y.value);
            result.derivativeOp = @ adr_divD;
            result.parents = [result.parents x y];
        end
        
        function modified_parents = adr_divD(prevDerivative, adNodes)
            adNodes(1).derivative = adNodes(1).derivative + prevDerivative*(1/adNodes(2).value);
            adNodes(2).derivative = adNodes(2).derivative + ...
                prevDerivative*(-adNodes(1).value/adNodes(2).value^2);
            modified_parents = adNodes;
        end
        
    end
   
end

















