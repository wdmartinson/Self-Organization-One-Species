function [IndicatorValue] = Indicator(distance,R)%arguments are distance and radius
eps=0.0001;
IndicatorValue =(tanh(((R-abs(distance))/eps))+1.)/2.;     % Smoothed indicator function   
end
