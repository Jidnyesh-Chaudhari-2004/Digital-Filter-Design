% Open a new Simulink model
new_system('system_diagram');

% Add blocks to the model
% Input block
in = add_block('system_diagram/Subsystem', 'simulink/Sources/Signal Generator');
set_param(in, 'Amplitude', '1', 'SampleTime', '1');

% Delay blocks
delay1 = add_block('system_diagram/Subsystem', 'simulink/Discrete/Delay');
set_param(delay1, 'Delay', '1');

delay2 = add_block('system_diagram/Subsystem', 'simulink/Discrete/Delay');
set_param(delay2, 'Delay', '2');

% Gain blocks
gain1 = add_block('system_diagram/Subsystem', 'simulink/Math Operations/Gain');
set_param(gain1, 'Gain', '-2*cos(theta)');

gain2 = add_block('system_diagram/Subsystem', 'simulink/Math Operations/Gain');
set_param(gain2, 'Gain', '1');

% Sum block
sum = add_block('system_diagram/Subsystem', 'simulink/Math Operations/Sum');
set_param(sum, 'Inputs', '3', 'IconShape', 'round');

% Output block
out = add_block('system_diagram/Subsystem', 'simulink/Sinks/Scope');

% Connect the blocks
add_line('system_diagram/Subsystem', in, delay1, 'autorouting','on');
add_line('system_diagram/Subsystem', in, gain2, 'autorouting','on');
add_line('system_diagram/Subsystem', in, sum, 'autorouting','on');
add_line('system_diagram/Subsystem', delay1, gain1, 'autorouting','on');
add_line('system_diagram/Subsystem', delay1, sum, 'autorouting','on');
add_line('system_diagram/Subsystem', delay2, sum, 'autorouting','on');
add_line('system_diagram/Subsystem', gain1, sum, 'autorouting','on');
add_line('system_diagram/Subsystem', gain2, sum, 'autorouting','on');
add_line('system_diagram/Subsystem', sum, out, 'autorouting','on');

% Open the model
open_system('system_diagram');