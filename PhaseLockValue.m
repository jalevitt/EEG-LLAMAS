function [PLV, Theta] = PhaseLockValue(SignalA, SignalB)
%Computes the phase lock value between two signals A and B
%  


N = length(SignalA);
HilA = hilbert(SignalA);
HilB = hilbert(SignalB);
PhaseA = angle(HilA);
PhaseB = angle(HilB);

Theta = abs(PhaseA - PhaseB);
Theta(Theta > pi) = 2 * pi - Theta(Theta > pi);
Theta = mean(Theta);
PLV = abs(sum(exp(1).^(1i * (PhaseA - PhaseB)))/N);
end

