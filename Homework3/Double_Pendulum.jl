using Plots

const m = 1.0
const g = 9.8
const l = 1.0

function dϕ1(ϕ1, ϕ2, p1, p2)
    1/(m*l^2) * (p1 - p2*cos(ϕ1 - ϕ2))/(1 + sin(ϕ1 - ϕ2)^2)
end

function dϕ2(ϕ1, ϕ2, p1, p2)
    1/(m*l^2) * (2*p2 - p1*cos(ϕ1 - ϕ2))/(1 + sin(ϕ1 - ϕ2)^2)
end

function dp1(ϕ1, ϕ2, p1, p2)
    1/(m*l^2) * 1/(1 + sin(ϕ1 - ϕ2)^2) * (-p1*p2*sin(ϕ1 - ϕ2) + (p1^2 + 2p2^2 - 2*p1*p2*cos(ϕ1 - ϕ2))/(1 + sin(ϕ1 - ϕ2)^2) * cos(ϕ1 - ϕ2) * sin(ϕ1 - ϕ2)) - 2*m*g*l*sin(ϕ1)
end

function dp2(ϕ1, ϕ2, p1, p2)
    1/(m*l^2) * 1/(1 + sin(ϕ1 - ϕ2)^2) * (p1*p2*sin(ϕ1 - ϕ2) - (p1^2 + 2p2^2 - 2*p1*p2*cos(ϕ1 - ϕ2))/(1 + sin(ϕ1 - ϕ2)^2) * cos(ϕ1 - ϕ2) * sin(ϕ1 - ϕ2)) - m*g*l*sin(ϕ2)
end

println("Done!")
