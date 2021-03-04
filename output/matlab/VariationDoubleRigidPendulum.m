function [varargout] = VariationDoubleRigidPendulum(varargin)
%% 
% This function was generated by the Scala: Dynamics on Manifold package
% Wed Mar 03 21:46:58 PST 2021
% 
% NOTE: This function makes use of geometric-toolbox
% 		 link: https://github.com/HybridRobotics/geometric-toolbox
%%


%% Function Inputs
params = varargin{1};

%% 
A = -[(((((((hat(l1))*(((-1.0))*(hat(((-1.0))*((transpose(R1))*((R2)*((hat(lc2))*((m2)*(((-1.0))*(dOmega_R2))))))))))+((m2)*((hat(l1))*(((-1.0))*(hat(((-1.0))*((transpose(R1))*((R2)*((hat(Omega_R2))*((hat(Omega_R2))*(lc2))))))))))))+(((((-1.0))*((((1.0))*((m1)*(g)))*((hat(lc1))*(((-1.0))*(hat(((-1.0))*((transpose(R1))*(e3))))))))+(((-1.0))*((((1.0))*((g)*(((mt)+(m2)))))*((hat(l1))*(((-1.0))*(hat(((-1.0))*((transpose(R1))*(e3))))))))))))+(((-1.0))*(((-1.0))*(hat(((-1.0))*((transpose(R1))*((R2)*(((-1.0))*(u))))))))),((((((hat(l1))*(transpose(R1)))*((R2)*(((-1.0))*(hat((hat(lc2))*((m2)*(((-1.0))*(dOmega_R2))))))))+((m2)*(((hat(l1))*(transpose(R1)))*((R2)*(((-1.0))*(hat((hat(Omega_R2))*((hat(Omega_R2))*(lc2))))))))))+(((-1.0))*((transpose(R1))*((R2)*(((-1.0))*(hat(((-1.0))*(u)))))))),(((hat(Omega_R1))*(((((J1)+(((-1.0))*((m1)*((hat(lc1))*(hat(lc1)))))))+(((-1.0))*((((m1)+(mt)))*((hat(l1))*(hat(l1))))))))+(((-1.0))*(hat((((((J1)+(((-1.0))*((m1)*((hat(lc1))*(hat(lc1)))))))+(((-1.0))*((((m1)+(mt)))*((hat(l1))*(hat(l1)))))))*(Omega_R1))))),(m2)*((((((hat(l1))*(transpose(R1)))*(R2))*(((-1.0))*(hat((hat(Omega_R2))*(lc2)))))+(((((hat(l1))*(transpose(R1)))*(R2))*(hat(Omega_R2)))*(((-1.0))*(hat(lc2)))))),;
((((hat(lc2))*(transpose(R2)))*((R1)*(((-1.0))*(hat((hat(l1))*((m2)*(((-1.0))*(dOmega_R1))))))))+((m2)*(((hat(lc2))*(transpose(R2)))*((R1)*(((-1.0))*(hat((hat(Omega_R1))*((hat(Omega_R1))*(l1))))))))),(((((hat(lc2))*(((-1.0))*(hat(((-1.0))*((transpose(R2))*((R1)*((hat(l1))*((m2)*(((-1.0))*(dOmega_R1))))))))))+((m2)*((hat(lc2))*(((-1.0))*(hat(((-1.0))*((transpose(R2))*((R1)*((hat(Omega_R1))*((hat(Omega_R1))*(l1))))))))))))+((((1.0))*((m2)*(g)))*((hat(lc2))*(((-1.0))*(hat(((-1.0))*((transpose(R2))*(e3)))))))),(m2)*((((((hat(lc2))*(transpose(R2)))*(R1))*(((-1.0))*(hat((hat(Omega_R1))*(l1)))))+(((((hat(lc2))*(transpose(R2)))*(R1))*(hat(Omega_R1)))*(((-1.0))*(hat(l1)))))),(((hat(Omega_R2))*(((J2)+(((-1.0))*((m2)*((hat(lc2))*(hat(lc2))))))))+(((-1.0))*(hat((((J2)+(((-1.0))*((m2)*((hat(lc2))*(hat(lc2)))))))*(Omega_R2))))),;
];

B = -[((-1.0))*(((-1.0))*((transpose(R1))*(R2))),;
((-1.0))*(eye(3)),;
];




%% Outputs
varargout{1} = A;
varargout{1} = B;

end