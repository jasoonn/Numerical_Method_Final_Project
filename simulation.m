clc;close all;clear all ;

dt = 0.00001; %every timestep stand for 0.01 second
speed = [0 -560 560];  %cm/s
ball_center = [0 675 230];
gravity_force = 980.6; % cm/s^2
mass = 0.625; % kg
area = 0.0454; %m^2
drag_coefficient = 0.49;
air_density = 1.2; % kg/m^3
magnus_cl = 0.15;
%building background
court_size = [600 900 600];
%line
for i =1:400
    base_line_lx(i) =-245;
    base_line_ly(i) = i;
    base_line_lz(i) = 0;
end
for i =1:580
    base_line_rx(i) =245;
    base_line_ry(i) = i;
    base_line_rz(i) = 0;
end
for i = -245:245
    base_line_midx(i+246) = i;
    base_line_midy(i+246)=580;
    base_line_midz(i+246)=0;
end
court_center = [0 580];
court_counter = 1;
for i =-245:245
    for j = 580:780
        if round(sqrt(i^2+(j-580)^2))==180
            base_roundx(court_counter) =i;
            base_roundy(court_counter) =j;
            base_roundz(court_counter) =0;
            court_counter= court_counter+1;
        end
    end
end
court_center = [0 157];
court_counter = 1;
for i =-300:300
    for j = 580:900
        if round(sqrt(i^2+(j-157.5)^2))==675  %675 is radius
            base_round_3x(court_counter) =i;
            base_round_3y(court_counter) =j;
            base_round_3z(court_counter) =0;
            court_counter= court_counter+1;
        end
    end
end

%rim
rim_counter = 1;
for i =-50:50
    for j = 0:157
        if abs(sqrt(i^2+(j-157.5)^2)-22.5)<0.1  %675 is radius
            rim_x(rim_counter) =i;
            rim_y(rim_counter) =j;
            rim_z(rim_counter) =305;
            rim_counter= rim_counter+1;
        end
    end
end
rim_counter = 1;
for i =-50:50
    for j = 158:300
        if abs(sqrt(i^2+(j-157.5)^2)-22.5)<0.1  %675 is radius
            rim_xu(rim_counter) =i;
            rim_yu(rim_counter) =j;
            rim_zu(rim_counter) =305;
            rim_counter= rim_counter+1;
        end
    end
end


%board
counter = 1;
for i =1:105
    board_x(counter) = -90;
    board_y(counter) = 120;
    board_z(counter) = 289+i;
    counter = counter+1;
    board_x(counter) = 90;
    board_y(counter) = 120;
    board_z(counter) = 289+i;
    counter = counter+1;
end


random_range = 0.04;

for decide_angle = 6:18
    for decide_omega = 1:16
        dt = 0.00001; 
        simulation_num = 4000;
        success = 0;
        fail = 0;
        rim_center = [0 157.5 305];
        Initial_angle = 35 + decide_angle*2;
        Initial_speed_magnitude = 1200;
        Changing = 36;
        last_sign = 1; %true: adding initial_speed    false: decreasing
        anguler_vector = [-0.2*pi-decide_omega*0.2 0 0];
        ball_radius = 12.3;

        %Finding Perfect trace under fix angle

        while true
            hitting = false;
            Initial_speed = [0 -Initial_speed_magnitude*cos(Initial_angle*pi/180) Initial_speed_magnitude*sin(Initial_angle*pi/180)];
            speed = Initial_speed;
            ball_center = [0 675 230];
            for i = 1:200000 % 0->5s
                ball_center = ball_center + speed*dt;
                ac = [0 0 -gravity_force];
                %drag_force
                drag_force = air_density * area * drag_coefficient * sum((speed/100).^2);
                drag_force = [0 -speed(2)/sqrt(sum(speed.^2)) -speed(3)/sqrt(sum(speed.^2))] .* drag_force;
                magnus_force = magnus_cl*16/3*pi^2*(ball_radius/100)^3*air_density* cross(anguler_vector,speed/100);
                ac = ac+ drag_force*100/mass + magnus_force*100/mass;
                speed = speed + ac*dt;
                if  mod(i,15000)==1  
                    a = ball(ball_center,12.3);
                    plot3(base_line_lx,-base_line_ly,base_line_lz,base_line_rx,-base_line_ry,base_line_rz,base_line_midx,-base_line_midy,base_line_midz,base_roundx,-base_roundy,base_roundz,base_round_3x,-base_round_3y,base_round_3z,rim_x,-rim_y,rim_z,rim_xu,-rim_yu,rim_zu,board_x,-board_y,board_z,a(1,:),-a(2,:),a(3,:))
                    axis equal
                    ylim([-1000 0])
                    zlim([0 600])
                    drawnow
                end
                if abs(ball_center(2)-rim_center(2))<0.01
                    if ball_center(3)-rim_center(3) >0.01
                        if last_sign == 1
                            Changing = Changing/2;
                        end
                        Initial_speed_magnitude = Initial_speed_magnitude - Changing;
                        last_sign = -1;
                    elseif ball_center(3)-rim_center(3) < -0.01
                        if last_sign == -1
                            Changing = Changing/2;
                        end
                        Initial_speed_magnitude = Initial_speed_magnitude + Changing;
                        last_sign = 1;
                    else 
                        hitting = true;
                    end
                    break;
                end

            end
            if hitting == true

                break;
            end
        end

        disp("Finished Finding, start simulation")


        for sim_num = 1:simulation_num
            dt = 0.00001; 
            disp(sim_num)
            rand_speed = rand*random_range*2 + (1-random_range);
            rand_angle = rand*random_range*2 + (1-random_range);
            rand_omega = rand*random_range*2 + (1-random_range);
            upper = false;%judging basketball
            ball_center = [0 675 230];
            ball_radius = 12.3;
            Initial_speed = [0 -Initial_speed_magnitude*cos(Initial_angle*rand_angle*pi/180) Initial_speed_magnitude*sin(Initial_angle*rand_angle*pi/180)]*rand_speed;
            Initial_omega = anguler_vector * rand_omega;
            speed = Initial_speed;
            rim_front = [0 180 305];
            rim_back = [0 135 305];

            counter = 0;
            while true
                if counter>600000
                    fail = fail + 1;
                    break;
                end
                if ball_center(3)>305
                    if upper ==false
                        upper = true;
                    end
                else
                    if upper == true
                        if 147<ball_center(2) && ball_center(2)<168
                            disp("hit")
                            success = success + 1;
                            disp(success/sim_num)
                        else
                            disp("miss")
                            disp(success/sim_num)
                            fail = fail + 1;
                        end
                        if mod(sim_num,100)==0
                            a = ball(ball_center,12.3);
                            plot3(base_line_lx,-base_line_ly,base_line_lz,base_line_rx,-base_line_ry,base_line_rz,base_line_midx,-base_line_midy,base_line_midz,base_roundx,-base_roundy,base_roundz,base_round_3x,-base_round_3y,base_round_3z,rim_x,-rim_y,rim_z,rim_xu,-rim_yu,rim_zu,board_x,-board_y,board_z,a(1,:),-a(2,:),a(3,:))
                            axis equal
                            ylim([-1000 0])
                            zlim([0 600])
                            drawnow
                        end
                        break
                    end
                end
                counter = counter + 1;
                ball_center = ball_center + speed*dt;
                ac = [0 0 -gravity_force];
                %drag_force
                drag_force = air_density * area * drag_coefficient * sum((speed/100).^2);
                drag_force = [0 -speed(2)/sqrt(sum(speed.^2)) -speed(3)/sqrt(sum(speed.^2))] .* drag_force;
                magnus_force = magnus_cl*16/3*pi^2*(ball_radius/100)^3*air_density* cross(Initial_omega,speed/100);
                %magnus_force = 0;

                ac = ac+ drag_force*100/mass + magnus_force*100/mass;%*100 because m/s^2 -> cm/s^2
                %if touch the board
                if ball_center(2) - ball_radius <= 120 && ball_center(3)>290  
                    dt = 0.000004;
                    normal_vector = [1 0];
                    force = touching_force(normal_vector, speed(2), 120 + ball_radius - ball_center(2))/mass;
                    ac = ac + [0 force(1) force(2)];
                end
                %if touch the front of the rim
                if sum((ball_center - rim_front).^2) <= ball_radius^2
                    dt = 0.000004;
                    normal_vector = (ball_center - rim_front)/ norm(ball_center - rim_front);
                    normal_speed = sum(normal_vector .* speed);
                    normal_vector = normal_vector(2:3);
                    force = touching_force(normal_vector, normal_speed, ball_radius - sqrt(sum((ball_center - rim_front).^2)))/mass;
                    
                    ac = ac + [0 force(1) force(2)];
                end
                %if touch the back of the rim
                %%{
                if sum((ball_center - rim_back).^2) <= ball_radius^2
                    dt = 0.000004;
                    normal_vector = (ball_center - rim_back)/ norm(ball_center - rim_back);
                    normal_speed = sum(normal_vector .* speed);
                    normal_vector = normal_vector(2:3);
                    force = touching_force(normal_vector, normal_speed, ball_radius - sqrt(sum((ball_center - rim_back).^2)))/mass;
                    ac = ac + [0 force(1) force(2)];
                end
                %}
                if  mod(counter,15000)==0 && mod(sim_num,100)==0
                    a = ball(ball_center,12.3);
                    plot3(base_line_lx,-base_line_ly,base_line_lz,base_line_rx,-base_line_ry,base_line_rz,base_line_midx,-base_line_midy,base_line_midz,base_roundx,-base_roundy,base_roundz,base_round_3x,-base_round_3y,base_round_3z,rim_x,-rim_y,rim_z,rim_xu,-rim_yu,rim_zu,board_x,-board_y,board_z,a(1,:),-a(2,:),a(3,:))
                    axis equal
                    ylim([-1000 0])
                    zlim([0 600])
                    %ylim([-300 0])
                    %zlim([100 400])
                    drawnow
                end

                speed = speed + ac*dt;    
            end
        end
        log(decide_angle,decide_omega) = success/simulation_num;
    end
end

function y = ball(center,radius)
    center = round(center);
    radius = round(radius);
    counter = 1;
    for i = center(1)-radius:center(1)+radius
        for j = center(2)-radius:center(2)+radius
            for k = center(3)-radius:center(3)+radius
                if abs(sqrt(sum(([i j k]-center).^2))-radius)<0.1
                    buff_x(counter) = i;
                    buff_y(counter) = j;
                    buff_z(counter) = k;
                    counter = counter + 1;
                end
            end
        end
    end
    y = [buff_x;buff_y;buff_z];     
end
function y = touching_force(normal_vector, normal_speed, normal_displacement)
k = 297600; %kg/cm
c = 227.664; %kg*s/cm
u = 0.1;
if c * normal_speed > 0.5 * normal_displacement * k 
    Nb = normal_displacement * k*0.5;
else
    Nb = normal_displacement * k - c * normal_speed;
end
Nb_vector = Nb * normal_vector;


Fb = u * Nb;
Fb = 0;
R = [cosd(90) -sind(90); sind(90) cosd(90)];


Fb_vector = Fb * (normal_vector*R);
y = Nb_vector + Fb_vector;
end
