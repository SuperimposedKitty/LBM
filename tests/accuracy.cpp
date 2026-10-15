#include "lbm/solver.hpp"
#include "lbm/solver3d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr double pi=3.14159265358979323846;
void check(bool ok,const char* message) {if (!ok) throw std::runtime_error(message);}

lbm::FoamSnapshot channel(int height, int width=8) {
    lbm::FoamSnapshot s;
    s.nx=width; s.ny=height+2; s.cells.resize(s.nx*s.ny);
    for (int y=0; y<s.ny; ++y)
        for (int x=0; x<s.nx; ++x) s.cells[y*s.nx+x].solid=y==0 || y==s.ny-1;
    return s;
}

double poiseuille(int h,bool interpolated) {
    auto s=channel(h);
    const double nu=0.1, lower=interpolated?0.25:0.5,
                 upper=h+(interpolated?0.35:0.5), height=upper-lower;
    const double umax=0.01*8/height, acceleration=8*nu*umax/(height*height);
    lbm::SolverConfig c;
    c.tau=0.8; c.flow.acceleration[0]=acceleration;
    c.flow.periodic_x=true; c.flow.interpolated_walls=interpolated;
    for (int y=0; y<s.ny; ++y)
        for (int x=0; x<s.nx; ++x) {
            auto& cell=s.cells[y*s.nx+x];
            if (interpolated) {
                cell.wall_distance=std::min(y-lower,upper-y);
                cell.wall_distance_known=true;
            }
            if (!cell.solid) cell.velocity[0]=acceleration*(y-lower)*(upper-y)/(2*nu);
        }
    lbm::Solver solver(s.nx,s.ny,c);
    solver.initialize_fields(s);
    const int steps=static_cast<int>(std::ceil(2*height*height/nu));
    const double initial_mass=solver.diagnostics().mass;
    double exchange=0;
    for (int t=0; t<steps; ++t) {
        solver.step_masked_flow(0);
        exchange+=solver.boundary_diagnostics().wall_mass_exchange;
    }
    double error=0,norm=0;
    for (int y=1; y<=h; ++y) {
        const double exact=acceleration*(y-lower)*(upper-y)/(2*nu);
        const double u=solver.grid().ux[solver.grid().scalar_index(3,y)];
        error+=(u-exact)*(u-exact); norm+=exact*exact;
    }
    check(std::abs(solver.diagnostics().mass-initial_mass-exchange)<1e-8,
          "Channel mass/operator balance failed");
    return std::sqrt(error/norm);
}

void channel_convergence() {
    for (bool interpolated:{false,true}) {
        double previous=0,last=0;
        for (int h:{8,16,32}) {
            const double error=poiseuille(h,interpolated);
            std::cout << "Poiseuille " << (interpolated?"BFL":"halfway") << " H=" << h
                      << " relativeL2=" << error;
            if (previous>0) {
                const double order=std::log(previous/error)/std::log(2.0);
                std::cout << " order=" << order;
                check(order>1.7,"Poiseuille observed order below 1.7");
            }
            std::cout << '\n'; previous=last=error;
        }
        check(last<0.01,"Fine channel velocity error exceeds 1 percent");
    }
}

void couette() {
    constexpr int h=16;
    auto s=channel(h);
    const double speed=0.02;
    for (int y=0; y<s.ny; ++y)
        for (int x=0; x<s.nx; ++x)
            s.cells[y*s.nx+x].velocity[0]=y==s.ny-1?speed:
                (y==0?0:speed*(y-0.5)/h);
    lbm::SolverConfig c; c.flow.periodic_x=true;
    lbm::Solver solver(s.nx,s.ny,c); solver.initialize_fields(s);
    for (int t=0; t<6000; ++t) solver.step_masked_flow(0);
    double error=0;
    for (int y=1; y<=h; ++y)
        error=std::max(error,std::abs(solver.grid().ux[solver.grid().scalar_index(3,y)]-
                                    speed*(y-0.5)/h)/speed);
    std::cout << "Couette max normalized error=" << error << '\n';
    check(error<1e-5,"Couette profile mismatch");
}

double taylor_green(int n) {
    lbm::FoamSnapshot s; s.nx=n; s.ny=n; s.cells.resize(n*n);
    const double amplitude=0.01*16/n, k=2*pi/n, nu=0.1;
    for (int y=0; y<n; ++y)
        for (int x=0; x<n; ++x) {
            auto& cell=s.cells[y*n+x];
            cell.velocity={-amplitude*std::cos(k*x)*std::sin(k*y),
                            amplitude*std::sin(k*x)*std::cos(k*y),0};
            cell.rho=1-0.75*amplitude*amplitude*(std::cos(2*k*x)+std::cos(2*k*y));
        }
    lbm::SolverConfig c; c.tau=0.8;
    lbm::Solver solver(n,n,c); solver.initialize_fields(s);
    const int steps=20*(n/16)*(n/16);
    for (int t=0; t<steps; ++t) solver.step();
    double error=0,norm=0;
    const double decay=std::exp(-2*nu*k*k*steps);
    for (int i=0; i<n*n; ++i)
        for (int d=0; d<2; ++d) {
            const double exact=s.cells[i].velocity[d]*decay;
            const double actual=d==0?solver.grid().ux[i]:solver.grid().uy[i];
            error+=(actual-exact)*(actual-exact); norm+=exact*exact;
        }
    return std::sqrt(error/norm);
}

void vortex_convergence() {
    double previous=0;
    for (int n:{16,32,64}) {
        const double error=taylor_green(n);
        std::cout << "TaylorGreen N=" << n << " relativeL2=" << error;
        if (previous>0) {
            const double order=std::log(previous/error)/std::log(2.0);
            std::cout << " order=" << order;
            check(order>1.7,"Taylor-Green observed order below 1.7");
        }
        if (n==64) check(error<0.01,"Fine Taylor-Green error exceeds 1 percent");
        std::cout << '\n'; previous=error;
    }
}

void forcing3d() {
    for (auto collision:{lbm::CollisionModel::BGK,lbm::CollisionModel::MRT}) {
        lbm::SolverConfig c; c.collision_model=collision;
        c.flow.force_density={2e-7,-1e-7,0}; c.flow.acceleration={1e-6,2e-6,0};
        lbm::FoamSnapshot s; s.nx=s.ny=6; s.cells.resize(36);
        for (auto& cell:s.cells) cell.rho=1.3;
        lbm::Solver solver(6,6,c); solver.initialize_fields(s); solver.run(20);
        check(std::abs(solver.grid().ux[0]-20*(2e-7/1.3+1e-6))<1e-12 &&
              std::abs(solver.grid().uy[0]-20*(-1e-7/1.3+2e-6))<1e-12,
              "2D Guo momentum integration mismatch");
        check(std::abs(solver.diagnostics().mass-1.3*36)<1e-10,"Forced 2D mass changed");
    }
    for (auto lattice:{lbm::Lattice3DModel::D3Q19,lbm::Lattice3DModel::D3Q27})
        for (auto collision:{lbm::CollisionModel::BGK,lbm::CollisionModel::MRT}) {
            lbm::Solver3DConfig c;
            c.lattice_model=lattice; c.collision_model=collision;
            c.flow.force_density={2e-7,-1e-7,3e-7};
            c.flow.acceleration={1e-6,2e-6,-1e-6};
            lbm::FoamSnapshot s; s.nx=s.ny=s.nz=6; s.two_dimensional=false; s.cells.resize(216);
            for (auto& cell:s.cells) cell.rho=1.3;
            lbm::Solver3D solver(6,6,6,c); solver.initialize_fields(s); solver.run(20);
            const auto& g=solver.grid();
            for (int d=0; d<3; ++d) {
                const double actual=d==0?g.ux[0]:(d==1?g.uy[0]:g.uz[0]);
                const double expected=20*(c.flow.force_density[d]/1.3+c.flow.acceleration[d]);
                check(std::abs(actual-expected)<1e-12,"3D Guo momentum integration mismatch");
            }
            check(std::abs(solver.diagnostics().mass-1.3*216)<1e-9,"Forced 3D mass changed");
        }
}

void boundary_moments() {
    for (const auto* l:{&lbm::detail::flow_lattice2d(),&lbm::d3q19,&lbm::d3q27}) {
        lbm::SinglePhaseOptions options;
        const double uz=l->q>9?0.001:0;
        options.force_density={1e-6,2e-6,l->q>9?1e-6:0};
        options.acceleration={2e-6,-3e-6,l->q>9?-2e-6:0};
        for (bool inlet:{true,false}) {
            std::array<double,27> f{};
            for (int q=0; q<l->q; ++q) f[q]=l->w[q]*(1+0.001*q);
            const auto before=f;
            lbm::detail::reconstruct_x_boundary(f,*l,options,inlet,{0.025,0.003,uz});
            double rho=0,jx=0,jy=0,jz=0;
            for (int q=0; q<l->q; ++q) {
                rho+=f[q]; jx+=l->cx[q]*f[q]; jy+=l->cy[q]*f[q]; jz+=l->cz[q]*f[q];
                if (l->cx[q]!=(inlet?1:-1)) check(f[q]==before[q],"Boundary changed a known population");
            }
            const auto force=lbm::detail::total_force(options,rho);
            check(std::abs((jy+0.5*force[1])/rho-0.003)<1e-12,"Boundary tangential velocity mismatch");
            check(std::abs((jz+0.5*force[2])/rho-uz)<1e-12,"Boundary z tangential velocity mismatch");
            if (inlet) check(std::abs((jx+0.5*force[0])/rho-0.025)<1e-12,"Velocity inlet mismatch");
            else check(std::abs(rho-1)<1e-12,"Pressure outlet mismatch");
        }
    }
}

double pipe_error(int radius,lbm::Lattice3DModel lattice=lbm::Lattice3DModel::D3Q19) {
    lbm::FoamSnapshot s;
    s.nx=4; s.ny=s.nz=2*radius+4; s.two_dimensional=false;
    s.cells.resize(s.nx*s.ny*s.nz);
    const double cy=0.5*(s.ny-1)+0.17, cz=0.5*(s.nz-1)-0.23;
    const double nu=0.1, speed=0.01*4/radius, acceleration=4*nu*speed/(radius*radius);
    for (int z=0; z<s.nz; ++z)
        for (int y=0; y<s.ny; ++y)
            for (int x=0; x<s.nx; ++x) {
                auto& cell=s.cells[(z*s.ny+y)*s.nx+x];
                const double r2=(y-cy)*(y-cy)+(z-cz)*(z-cz);
                cell.solid=r2>=radius*radius;
                cell.wall_distance=radius-std::sqrt(r2); cell.wall_distance_known=true;
                if (!cell.solid) cell.velocity[0]=acceleration*(radius*radius-r2)/(4*nu);
            }
    lbm::Solver3DConfig c;
    c.lattice_model=lattice; c.flow.periodic_x=true; c.flow.interpolated_walls=true;
    c.flow.acceleration[0]=acceleration;
    lbm::Solver3D solver(s.nx,s.ny,s.nz,c); solver.initialize_fields(s);
    const int steps=static_cast<int>(2*radius*radius/nu);
    for (int t=0; t<steps; ++t) solver.step_masked_flow(0);
    double error=0,norm=0;
    for (std::size_t i=0; i<s.cells.size(); ++i) if (!s.cells[i].solid) {
        const double exact=s.cells[i].velocity[0], diff=solver.grid().ux[i]-exact;
        error+=diff*diff; norm+=exact*exact;
    }
    return std::sqrt(error/norm);
}

void curved_pipe_convergence() {
    double last=0;
    for (int radius:{4,8,16}) {
        const double error=pipe_error(radius);
        std::cout << "3D circular pipe R=" << radius << " relativeL2=" << error;
        if (last>0) std::cout << " order=" << std::log(last/error)/std::log(2.0);
        std::cout << '\n';
        if (last>0) check(error<last,"Curved pipe error did not decrease");
        last=error;
    }
    check(last<0.005,"Curved pipe fine velocity error exceeds 0.5 percent");
    check(pipe_error(8,lbm::Lattice3DModel::D3Q27)<0.03,"D3Q27 curved pipe mismatch");
}

void cavity_reference() {
    // Ghia Re=100 中心速度，见 Sandia Fuego Verification Manual 第 22 章。
    for (int n:{32,48}) {
        const double speed=0.04*32/n, nu=speed*n/100;
        lbm::SolverConfig c; c.tau=0.5+3*nu;
        lbm::Solver solver(n+2,n+2,c); solver.initialize_lid_driven_cavity(speed);
        const int steps=static_cast<int>(40*n/speed);
        solver.run_lid_driven_cavity(steps,speed);
        const auto& g=solver.grid();
        double u=0,v=0;
        for (int y=n/2; y<=n/2+1; ++y)
            for (int x=n/2; x<=n/2+1; ++x) {
                const auto i=g.scalar_index(x,y); u+=g.ux[i]/(4*speed); v+=g.uy[i]/(4*speed);
            }
        const double error=std::max(std::abs(u+0.20581),std::abs(v-0.05454));
        std::cout << "Cavity Re=100 N=" << n << " center U,V=" << u << ',' << v
                  << " max absolute normalized error=" << error << '\n';
        check(error<0.01,"Cavity center velocity disagrees with Ghia reference");
    }
}

void cylinder_reference(int diameter) {
    check(diameter>=10 && diameter<=40 && diameter%10==0,"Cylinder diameter must be 10,20,30,40");
    // DFG 2D-1 几何与抛物入口；本实现出口采用固定参考密度，参考为 do-nothing 出口。
    const int nx=22*diameter, h=41*diameter/10;
    auto s=channel(h,nx);
    const double cx=2.0*diameter,cy=2.0*diameter,radius=0.5*diameter;
    const double mean_speed=0.02*10/diameter,nu=mean_speed*diameter/20;
    for (int y=0; y<s.ny; ++y)
        for (int x=0; x<s.nx; ++x) {
            auto& cell=s.cells[y*s.nx+x];
            const double xx=x+0.5,yy=y-0.5;
            cell.wall_distance=std::hypot(xx-cx,yy-cy)-radius; cell.wall_distance_known=true;
            cell.obstacle=cell.wall_distance<=0; cell.solid=cell.solid || cell.obstacle;
            const double eta=yy/h;
            cell.inlet_scale=std::max(0.0,6*eta*(1-eta));
            if (!cell.solid) {
                cell.velocity[0]=mean_speed*cell.inlet_scale;
                cell.boundary=x==0?1:(x==nx-1?2:0);
            }
        }
    lbm::SolverConfig c; c.tau=0.5+3*nu;
    c.flow.open_boundary=lbm::OpenBoundaryModel::VelocityPressure; c.flow.interpolated_walls=true;
    lbm::Solver solver(s.nx,s.ny,c); solver.initialize_fields(s);
    const int steps=30000*(diameter/10)*(diameter/10);
    double old_drag=0,drag=0,lift=0;
    for (int t=1; t<=steps; ++t) {
        solver.step_masked_flow(mean_speed);
        if (t%10000==0 || t==steps) {
            const auto& force=solver.boundary_diagnostics().obstacle_force;
            drag=2*force[0]/(mean_speed*mean_speed*diameter);
            lift=2*force[1]/(mean_speed*mean_speed*diameter);
            std::cout << "Cylinder Re=20 D=" << diameter << " step=" << t
                      << " Cd=" << drag << " Cl=" << lift << " deltaCd=" << drag-old_drag << '\n';
            old_drag=drag;
        }
    }
    const double error=std::abs(drag/5.57953523384-1);
    std::cout << "Cylinder reference Cd=5.57953523384 relative error=" << error << '\n';
    check(error<0.05,"Cylinder drag error exceeds 5 percent; see resolution/outlet differences");
    check(std::abs(lift-0.010618948146)<0.02,"Cylinder lift differs excessively from reference");
}
} // namespace

int main(int argc,char** argv) {
    try {
        std::cout << std::unitbuf;
        if (argc>1) {
            const bool reference=std::string(argv[1])=="--reference";
            check((reference || std::string(argv[1])=="--cylinder") && argc<=3,
                  "Usage: lbm_accuracy_tests [--reference|--cylinder [diameter]]");
            int diameter=10;
            if (argc==3) {std::size_t used=0; diameter=std::stoi(argv[2],&used); check(used==std::string(argv[2]).size(),"Invalid diameter");}
            if (reference) {curved_pipe_convergence(); cavity_reference();}
            cylinder_reference(diameter);
        } else {
            boundary_moments(); forcing3d(); couette(); channel_convergence(); vortex_convergence();
        }
        std::cout << "Single-phase accuracy benchmarks passed\n";
        return 0;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
