

data = expData(1);

t = data.t;
u = data.u;
theta1 = data.theta_1;
theta2 = data.theta_2;

simin = [t, u];

theta_1_val = [t, theta1];
theta_2_val = [t, theta2];

theta_1_0_compare = theta1(1);
theta_2_0_compare = theta2(1);