using GLMakie
GLMakie.activate!(inline=false)

function x_deriv(y,t) # Submit as [x', v_x']
    v_p = y[1]
    v = y[2]
    g = 9.81
    return [ v , 0 ]
end

function y_deriv(y,t) # Submit as [y', v_y']
    v_p = y[1]
    v = y[2]
    g = 9.81
    return [ v , 0 ]
end

y0_x = [2,2]
y0_y = [2,10]

#st_time = 0
#end_time = 5
#h = 10^-2 #step size
#r = (end_time - st_time)/h
#t = LinRange(st_time, end_time, trunc(Int,r)+1)
t0 = 2.3
function rk4(f,y0,t0) #t0 initial time, h-time step
    args = 1
    r = 20
    boundary_const = 200
    h = 10^(-3)
    t = [t0, t0+h]
    N = size(t)[1]
    y = zeros(N,2)
    y[1,:] = y0
    for i=1:N-1
        h = t[i+1] - t[i]
        k1 = f(y[i,:],t[i])
        k2 = f(y[i,:] + k1*h/2, t[i] + h/2)
        k3 = f(y[i,:] + k2*h/2, t[i] + h/2)
        k4 = f(y[i,:] + k3*h, t[i] + h)
        y[i+1,:] = y[i,:] + h*(k1+2*k2+2*k3+k4)/6
    end
    if y[2,1] >=  boundary_const || y[2,1] <= 0
        y[2,2] = -y[2,2]
    end
    return y[2,:]
end


sol_x = rk4(x_deriv, y0_x, t0)
sol_y = rk4(y_deriv, y0_y, t0)

x_dat = vec(sol_x[1:end,[1]])
v_x_dat = vec(sol_x[1:end,[2]]) 

y_dat = vec(sol_y[1:end,[1]])
v_y_dat = vec(sol_y[1:end,[2]])

plot(x_dat,y_dat)

# Circle(Point2f(0, 0), 15f0) ####Circle plotting in GLMakie

# Point2f
# rand(Point2f, 10) 


mutable struct Hard_Ball
    Position::Point2f
    Velocity::Vector{Float64}
    radius::Float64         
end


# typeof(Point2f(1,1))


function norm(p,q) #distance function
    return sqrt(sum((p .- q).^2))
end

# getindex()

# for q in particles
#     if norm(p,q) < r
#         push!(neighbour_arr,q)
#         push!(neighbour_vel_arr, vels[getindex(q)])
#     end 
# end

function neighbour_find(r,p,particles,vels) #neighbour of p in particle at radius r 
    neighbour_arr = Point{2, Float32}[]
    neighbour_vel_arr = Point{2, Float32}[]
    for i = 1:size(particles)[1]
        if norm(p,particles[i]) <= r
            push!(neighbour_arr,particles[i]) 
            push!(neighbour_vel_arr, vels[i])
        end
    end
    return neighbour_arr,neighbour_vel_arr
end




# particle_1 = Hard_Ball([1.2,2.3],[1,1],0.2)

# list_parts = [Hard_Ball(rand(Point2f),rand(Float64,2),rand(Float64)) for i = 1:10]

# [i.Position for i in list_parts]

# list_parts[1].Position

# rand(Point2f)
# rand(0:0.01:0.3,1)
# rand(Float64,2)

function dot(x,y)
    return sum(x.*y)
end

#dot([1,2,1],[1,2,3])

function collide_vel(p, q, vp, vq) #p and q , vp and vq .... the return velocities are vp' and vq' reply
    return @. p-0.07, q+0.07, vp - (dot(vp-vq,p-q)/(norm(p,q)^2)) * (p-q) , vq - (dot(vq-vp,q-p)/(norm(q,p)^2)) * (q-p)
end


function particles_data_update(part_pos,part_vel,part_rad,t) #part_pos = particle_position
    #rad = 0.02 #radius of the discs ### Old
    for i=1:size(part_pos)[1]
        for j=1:size(part_pos)[1]
            if norm(part_pos[j],part_pos[i]) <= part_rad[i] + part_rad[j] && j!=i
                (part_vel[i], part_vel[j]) = collide_vel(part_pos[i], part_pos[j], part_vel[i], part_vel[j])
                (part_pos[i],part_pos[j]) = (part_pos[i] .+0.01,part_pos[j] .-0.01) ###########
                # diff = part_pos[i] .- part_pos[j]
                # (part_pos[i],part_pos[j]) = (part_pos[i] .- diff ,part_pos[j] .+ diff) ###########
            end
        end
    end
    #println(part_vel)

    part_pos = [Point2f(rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[1],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[1] ) for i=1:size(part_pos)[1]]
    part_vel = [[rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[2],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[2] ] for i=1:size(part_pos)[1]]

    return part_pos, part_vel

end

function particles_data_update_new(part_pos,part_vel,part_rad,t) #part_pos = particle_position
    #rad = 0.02 #radius of the discs ### Old
        for i=1:size(part_pos)[1]
            for j=1:size(part_pos)[1]
                if norm(part_pos[j],part_pos[i]) <= 36 && j!=i #part_rad[i] + part_rad[j]
                (part_pos[i],part_pos[j],part_vel[i], part_vel[j]) = collide_vel(part_pos[i], part_pos[j], part_vel[i], part_vel[j])
                # (part_pos[i],part_pos[j]) = (part_pos[i] .+0.01,part_pos[j] .-0.01) ###########
                # diff = part_pos[i] .- part_pos[j]
                # (part_pos[i],part_pos[j]) = (part_pos[i] .- diff ,part_pos[j] .+ diff) ###########
                end
            end
        # part_pos = [Point2f(rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[1],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[1] ) for i=1:size(part_pos)[1]]
        # part_vel = [[rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[2],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[2] ] for i=1:size(part_pos)[1]]    
    end
    #println(part_vel)

    part_pos = [Point2f(rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[1],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[1] ) for i=1:size(part_pos)[1]]
    part_vel = [[rk4(x_deriv, [part_pos[i][1], part_vel[i][1]], t)[2],rk4(y_deriv, [part_pos[i][2], part_vel[i][2]], t)[2] ] for i=1:size(part_pos)[1]]
    
    return part_pos, part_vel
end

list_particles = [Hard_Ball(rand(Point2f).* Point2f(200, 200),rand(Float64,2).*400,20) for i = 1:3]

function main()
        t_sys = 0
        s = Scene(resolution = (200, 200), camera = campixel!)
        list_particles = [Hard_Ball(rand(Point2f).* Point2f(200, 200),rand(Float64,2).*500,18) for i = 1:3]

        pos_list = Observable([p.Position for p in list_particles[1:end]])
        vel_list = Observable([p.Velocity for p in list_particles[1:end]])
        rad_list = [p.radius for p in list_particles[1:end]]


        # data = Observable(rand(Point2f, 10_000) .* Point2f(800, 800))
        scatter!(s, pos_list,markersize=40 )#[p.Position for p in list_particles[][1:end]])
        screen = display(s)
        task = @async while isopen(screen)
            # data[] .+= randn.(Point2f) # update rule
            Hehehe = particles_data_update_new(pos_list[], vel_list[], rad_list, t_sys )
    
            for i in 1:size(list_particles)[1]
                # list_particles[i].Position = Hehehe[1][i]
                # list_particles[i].Velocity = Hehehe[2][i]
                pos_list[][i] = Hehehe[1][i]
                vel_list[][i] = Hehehe[2][i]
            end
            # notify(list_particles)
            notify(pos_list), notify(vel_list)
            sleep(1/100000)
            t_sys+=0.001
            #print(t_sys)
            #notify(t_sys)
        end
        Base.errormonitor(task)        
end

main()


function main_new()
    t_sys = 0
    #s = Scene(resolution = (400, 400), camera = campixel!)
    f = Figure()
    Axis(f[1, 1], aspect = DataAspect())
    list_particles = [Hard_Ball(rand(Point2f).* Point2f(400, 400),rand(Float64,2).*30,5) for i = 1:10]

    pos_list = Observable([p.Position for p in list_particles[1:end]])
    vel_list = Observable([p.Velocity for p in list_particles[1:end]])
    rad_list = [p.radius for p in list_particles[1:end]]

    for i=1:size(list_particles)[1]
        poly!(Circle(pos_list[][i] , 5), color = :pink)
    end
    # data = Observable(rand(Point2f, 10_000) .* Point2f(800, 800))
    #scatter!(s, pos_list)#[p.Position for p in list_particles[][1:end]])
    screen = display(f)
    task = @async while isopen(screen)
        # data[] .+= randn.(Point2f) # update rule
        Hehehe = particles_data_update(pos_list[], vel_list[], rad_list, t_sys )

        for i in 1:size(list_particles)[1]
            # list_particles[i].Position = Hehehe[1][i]
            # list_particles[i].Velocity = Hehehe[2][i]
            pos_list[][i] = Hehehe[1][i]
            vel_list[][i] = Hehehe[2][i]
        end
        # notify(list_particles)
        notify(pos_list), notify(vel_list)
        for i=1:size(list_particles)[1]
            poly!(Circle(pos_list[][i] , 5), color = :pink)
        end
        sleep(1/100000)
        t_sys+=0.01
        print(t_sys)
        #notify(t_sys)
    end
    Base.errormonitor(task)        
end


main_new()

# begin
#     s = Scene(resolution = (800, 800), camera = campixel!)
#     data = Observable(rand(Point2f, 10_000) .* Point2f(800, 800))
#     scatter!(s, data, color = rand(10_000))
#     screen = display(s)
#     task = @async while isopen(screen)
#         data[] .+= randn.(Point2f)
#         notify(data)
#         sleep(1/100)
#     end
#     Base.errormonitor(task)
# end|