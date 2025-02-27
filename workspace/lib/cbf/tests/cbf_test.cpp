#include <iostream>
#include <ginac/ginac.h>

#define M_PI 3.14159265358979323846264338327950288

int main()
{
    int STATE_VARS = 6;
    int INPUT_VARS = 3;
    
    GiNaC::symbol px("px"), py("py"), thx("thx"), vx("vx"), vy("vy"), w("w"), xt("xt"), yt("yt"), Ds("Ds"), fov("fov"), Rs("Rs");
    GiNaC::ex poly;


    GiNaC::matrix A = {{0, 0, 0, 1, 0, 0},
                    {0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 1},
                    {0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0}};

    GiNaC::matrix B = {{0, 0, 0},
                    {0, 0, 0},
                    {0, 0, 0},
                    {1, 0, 0},
                    {0, 1, 0},
                    {0, 0, 1}};

    
    GiNaC::matrix x = {{px, py, thx, vx, vy, w}};
    x = x.transpose();
    GiNaC::matrix x_target = {{xt, yt}};
    x_target = x_target.transpose();

    GiNaC::matrix f = A.mul(x);
    GiNaC::matrix g = B;
    std::cout << "f: " << f << std::endl;
    
    /* ----------- SAFETY CBF ---------- */
    GiNaC::ex norm2 = GiNaC::pow(x_target[0]-x[0], 2) + GiNaC::pow(x_target[1]-x[1], 2);
    GiNaC::ex b1 = norm2 - GiNaC::pow(Ds, 2);
    std::cout << "b1 : " << b1 << std::endl;

    GiNaC::matrix grad_b1 = GiNaC::matrix(STATE_VARS, 1);
    grad_b1(0, 0) = GiNaC::diff(b1, px);
    grad_b1(1, 0) = GiNaC::diff(b1, py);
    grad_b1(2, 0) = GiNaC::diff(b1, thx);
    grad_b1(3, 0) = GiNaC::diff(b1, vx);
    grad_b1(4, 0) = GiNaC::diff(b1, vy);
    grad_b1(5, 0) = GiNaC::diff(b1, w);
    std::cout << "Grad b: " << grad_b1 << std::endl;
    std::cout << "gradb shape: " << grad_b1.rows() << grad_b1.cols() << std::endl;
    std::cout << "f shape: " << f.rows() << f.cols() << std::endl;

    GiNaC::ex lfb1 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lfb1 = lfb1 + grad_b1(i, 0) * f(i, 0);
    }
    std::cout << "Lfb1: " << lfb1 << std::endl;
    
    GiNaC::matrix grad2_b1 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b1(0, 0) = GiNaC::diff(grad_b1, px); 
    grad2_b1(1, 0) = GiNaC::diff(grad_b1, py); 
    grad2_b1(2, 0) = GiNaC::diff(grad_b1, thx); 
    grad2_b1(3, 0) = GiNaC::diff(grad_b1, vx); 
    grad2_b1(4, 0) = GiNaC::diff(grad_b1, vy); 
    grad2_b1(5, 0) = GiNaC::diff(grad_b1, w);
    GiNaC::ex lf2b1 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lf2b1 = lf2b1 + grad2_b1(i, 0) * f(i, 0);
    }
    std::cout << "Lf2b1: " << lf2b1 << std::endl;

    GiNaC::matrix Ac1 = {{2*px - 2*xt, 2*py - 2*yt, 0}};
    GiNaC::ex B1 = lf2b1;
    GiNaC::ex B2 = 2 * b1 * lfb1;
    GiNaC::ex B3 = GiNaC::pow(lfb1, 2);
    GiNaC::ex B4 = 2 * GiNaC::pow(b1, 2) * lfb1;
    GiNaC::ex B5 = GiNaC::pow(b1, 4);

    GiNaC::ex Bc1 = B1 + B2 + B3 + B4 + B5;
    
    std::cout << "Ac1: " << Ac1 << std::endl;
    std::cout << "Bc1: " << Bc1.evalm() << std::endl;
    std::cout << std::endl << std::endl;

    /*
    GiNaC::ex eval_bc = GiNaC::subs(Bc1, px==0.0);
    eval_bc = GiNaC::subs(eval_bc, py==0.0);
    eval_bc = GiNaC::subs(eval_bc, thx==GiNaC::Pi/4);
    eval_bc = GiNaC::subs(eval_bc, vx==0.0);
    eval_bc = GiNaC::subs(eval_bc, vy==0.0);
    eval_bc = GiNaC::subs(eval_bc, w==0.0);
    eval_bc = GiNaC::subs(eval_bc, Ds==2.0);
    eval_bc = GiNaC::subs(eval_bc, xt==4.0);
    eval_bc = GiNaC::subs(eval_bc, yt==4.0);
    std::cout << "Bc: " << eval_bc << std::endl;

    GiNaC::ex eval_a = GiNaC::subs(Ac1, px==0.0);
    eval_a = GiNaC::subs(eval_a, py==0.0);
    eval_a = GiNaC::subs(eval_a, xt==4.0);
    eval_a = GiNaC::subs(eval_a, yt==4.0);
    std::cout << "Ac: " << eval_a << std::endl;
    */

    /* ----------- BORDER 1 CBF ---------- */
    GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});
    GiNaC::matrix R = {{GiNaC::cos(thx), GiNaC::sin(thx)},
                        {-GiNaC::sin(thx), GiNaC::cos(thx)}};
    GiNaC::matrix xt_rel = R.mul(d2);

    GiNaC::ex b2 = GiNaC::tan(fov/2)*xt_rel(0,0) + xt_rel(1,0);
    GiNaC::matrix grad_b2 = GiNaC::matrix(STATE_VARS, 1);
    grad_b2(0, 0) = GiNaC::diff(b2, px);
    grad_b2(1, 0) = GiNaC::diff(b2, py);
    grad_b2(2, 0) = GiNaC::diff(b2, thx);
    grad_b2(3, 0) = GiNaC::diff(b2, vx);
    grad_b2(4, 0) = GiNaC::diff(b2, vy);
    grad_b2(5, 0) = GiNaC::diff(b2, w);
    std::cout << "Grad b2: " << grad_b2 << std::endl;
    std::cout << "gradb shape: " << grad_b2.rows() << grad_b2.cols() << std::endl;
    std::cout << "f shape: " << f.rows() << f.cols() << std::endl;

    GiNaC::ex lfb2 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lfb2 = lfb2 + grad_b2(i, 0) * f(i, 0);
    }
    std::cout << "Lfb2: " << lfb2 << std::endl;

    GiNaC::matrix grad2_b2 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b2(0, 0) = GiNaC::diff(grad_b2, px); 
    grad2_b2(1, 0) = GiNaC::diff(grad_b2, py); 
    grad2_b2(2, 0) = GiNaC::diff(grad_b2, thx); 
    grad2_b2(3, 0) = GiNaC::diff(grad_b2, vx); 
    grad2_b2(4, 0) = GiNaC::diff(grad_b2, vy); 
    grad2_b2(5, 0) = GiNaC::diff(grad_b2, w);
    GiNaC::ex lf2b2 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lf2b2 = lf2b2 + grad2_b2(i, 0) * f(i, 0);
    }
    std::cout << "Lf2b2: " << lf2b2 << std::endl;

    GiNaC::matrix Ac2 = {{GiNaC::tan(fov/2), 1, -xt_rel(1,0)*GiNaC::tan(fov/2)+xt_rel(0,0)}};
    GiNaC::ex B12 = lf2b2;
    GiNaC::ex B22 = 2 * b2 * lfb2;
    GiNaC::ex B32 = GiNaC::pow(lfb2, 2);
    GiNaC::ex B42 = 2 * GiNaC::pow(b2, 2) * lfb2;
    GiNaC::ex B52 = GiNaC::pow(b2, 4);

    GiNaC::ex Bc2 = B12 + B22 + B32 + B42 + B52;
    
    std::cout << "Ac2: " << Ac2 << std::endl;
    std::cout << "Bc2: " << Bc2.evalm() << std::endl;
    std::cout << std::endl << std::endl;


    // TEST
    /*
    GiNaC::ex beta = 120.0 * GiNaC::Pi / 180.0;
    GiNaC::ex eval_bc = GiNaC::subs(Bc2, px==0.0);
    eval_bc = GiNaC::subs(eval_bc, py==0.0);
    eval_bc = GiNaC::subs(eval_bc, thx==GiNaC::Pi/4);
    eval_bc = GiNaC::subs(eval_bc, vx==0.0);
    eval_bc = GiNaC::subs(eval_bc, vy==0.0);
    eval_bc = GiNaC::subs(eval_bc, w==0.0);
    eval_bc = GiNaC::subs(eval_bc, Ds==2.0);
    eval_bc = GiNaC::subs(eval_bc, xt==4.0);
    eval_bc = GiNaC::subs(eval_bc, yt==4.0);
    eval_bc = GiNaC::subs(eval_bc, fov==beta);
    std::cout << "Bc: " << eval_bc.evalf() << std::endl;

    GiNaC::ex eval_a = GiNaC::subs(Ac2, px==0.0);
    eval_a = GiNaC::subs(eval_a, py==0.0);
    eval_a = GiNaC::subs(eval_a, xt==4.0);
    eval_a = GiNaC::subs(eval_a, yt==4.0);
    eval_a = GiNaC::subs(eval_a, fov==beta);
    eval_a = GiNaC::subs(eval_a, thx==GiNaC::Pi/4);
    std::cout << "Ac: " << eval_a.evalf() << std::endl;
    */


    /* ------------ BORDER 2 CBF ----------- */
    GiNaC::ex b3 = GiNaC::tan(fov/2)*xt_rel(0,0) - xt_rel(1,0);
    GiNaC::matrix grad_b3 = GiNaC::matrix(STATE_VARS, 1);
    grad_b3(0, 0) = GiNaC::diff(b3, px);
    grad_b3(1, 0) = GiNaC::diff(b3, py);
    grad_b3(2, 0) = GiNaC::diff(b3, thx);
    grad_b3(3, 0) = GiNaC::diff(b3, vx);
    grad_b3(4, 0) = GiNaC::diff(b3, vy);
    grad_b3(5, 0) = GiNaC::diff(b3, w);
    std::cout << "Grad b3: " << grad_b3 << std::endl;
    std::cout << "gradb shape: " << grad_b3.rows() << grad_b3.cols() << std::endl;
    std::cout << "f shape: " << f.rows() << f.cols() << std::endl;

    GiNaC::ex lfb3 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lfb3 = lfb3 + grad_b3(i, 0) * f(i, 0);
    }
    std::cout << "Lfb3: " << lfb3 << std::endl;

    GiNaC::matrix grad2_b3 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b3(0, 0) = GiNaC::diff(grad_b3, px); 
    grad2_b3(1, 0) = GiNaC::diff(grad_b3, py); 
    grad2_b3(2, 0) = GiNaC::diff(grad_b3, thx); 
    grad2_b3(3, 0) = GiNaC::diff(grad_b3, vx); 
    grad2_b3(4, 0) = GiNaC::diff(grad_b3, vy); 
    grad2_b3(5, 0) = GiNaC::diff(grad_b3, w);
    GiNaC::ex lf2b3 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lf2b3 = lf2b3 + grad2_b3(i, 0) * f(i, 0);
    }
    std::cout << "Lf2b3: " << lf2b3 << std::endl;

    GiNaC::matrix Ac3 = {{GiNaC::tan(fov/2), 1, -xt_rel(1,0)*GiNaC::tan(fov/2)-xt_rel(0,0)}};
    GiNaC::ex B13 = lf2b3;
    GiNaC::ex B23 = 2 * b3 * lfb3;
    GiNaC::ex B33 = GiNaC::pow(lfb3, 2);
    GiNaC::ex B43 = 2 * GiNaC::pow(b3, 2) * lfb3;
    GiNaC::ex B53 = GiNaC::pow(b3, 4);

    GiNaC::ex Bc3 = B13 + B23 + B33 + B43 + B53;
    
    std::cout << "Ac3: " << Ac3 << std::endl;
    std::cout << "Bc3: " << Bc3.evalm() << std::endl;
    std::cout << std::endl << std::endl;

    /* ----------- MAX DISTANCE CBF ------------ */
    GiNaC::ex b4 = -norm2 + GiNaC::pow(Rs, 2);
    std::cout << "b4 : " << b4 << std::endl;

    GiNaC::matrix grad_b4 = GiNaC::matrix(STATE_VARS, 1);
    grad_b4(0, 0) = GiNaC::diff(b4, px);
    grad_b4(1, 0) = GiNaC::diff(b4, py);
    grad_b4(2, 0) = GiNaC::diff(b4, thx);
    grad_b4(3, 0) = GiNaC::diff(b4, vx);
    grad_b4(4, 0) = GiNaC::diff(b4, vy);
    grad_b4(5, 0) = GiNaC::diff(b4, w);
    std::cout << "Grad b: " << grad_b4 << std::endl;
    std::cout << "gradb shape: " << grad_b4.rows() << grad_b4.cols() << std::endl;
    std::cout << "f shape: " << f.rows() << f.cols() << std::endl;

    GiNaC::ex lfb4 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lfb4 = lfb4 + grad_b4(i, 0) * f(i, 0);
    }
    std::cout << "Lfb4: " << lfb4 << std::endl;
    
    GiNaC::matrix grad2_b4 = GiNaC::matrix(STATE_VARS, 1);
    grad2_b4(0, 0) = GiNaC::diff(grad_b4, px); 
    grad2_b4(1, 0) = GiNaC::diff(grad_b4, py); 
    grad2_b4(2, 0) = GiNaC::diff(grad_b4, thx); 
    grad2_b4(3, 0) = GiNaC::diff(grad_b4, vx); 
    grad2_b4(4, 0) = GiNaC::diff(grad_b4, vy); 
    grad2_b4(5, 0) = GiNaC::diff(grad_b4, w);
    GiNaC::ex lf2b4 = 0.0;
    for (int i = 0; i < STATE_VARS; i++)
    {
        lf2b4 = lf2b4 + grad2_b4(i, 0) * f(i, 0);
    }
    std::cout << "Lf2b4: " << lf2b4 << std::endl;

    GiNaC::matrix Ac4 = {{-2*px + 2*xt, -2*py + 2*yt, 0}};
    GiNaC::ex B14 = lf2b4;
    GiNaC::ex B24 = 2 * b4 * lfb4;
    GiNaC::ex B34 = GiNaC::pow(lfb4, 2);
    GiNaC::ex B44 = 2 * GiNaC::pow(b4, 2) * lfb4;
    GiNaC::ex B54 = GiNaC::pow(b4, 4);

    GiNaC::ex Bc4 = B14 + B24 + B34 + B44 + B54;
    
    std::cout << "Ac4: " << Ac4 << std::endl;
    std::cout << "Bc4: " << Bc4.evalm() << std::endl;
    std::cout << std::endl << std::endl;

    /*-------------- TEST -----------
    GiNaC::ex beta = 120.0 * GiNaC::Pi / 180.0;
    GiNaC::ex eval_bc = GiNaC::subs(Bc4, px==0.0);
    eval_bc = GiNaC::subs(eval_bc, py==0.0);
    eval_bc = GiNaC::subs(eval_bc, thx==GiNaC::Pi/4);
    eval_bc = GiNaC::subs(eval_bc, vx==0.0);
    eval_bc = GiNaC::subs(eval_bc, vy==0.0);
    eval_bc = GiNaC::subs(eval_bc, w==0.0);
    eval_bc = GiNaC::subs(eval_bc, Ds==2.0);
    eval_bc = GiNaC::subs(eval_bc, xt==4.0);
    eval_bc = GiNaC::subs(eval_bc, yt==4.0);
    eval_bc = GiNaC::subs(eval_bc, fov==beta);
    eval_bc = GiNaC::subs(eval_bc, Rs==10.0);
    std::cout << "Bc: " << eval_bc.evalf() << std::endl;

    GiNaC::ex eval_a = GiNaC::subs(Ac4, px==0.0);
    eval_a = GiNaC::subs(eval_a, py==0.0);
    eval_a = GiNaC::subs(eval_a, xt==4.0);
    eval_a = GiNaC::subs(eval_a, yt==4.0);
    eval_a = GiNaC::subs(eval_a, fov==beta);
    eval_a = GiNaC::subs(eval_a, thx==GiNaC::Pi/4);
    eval_a = GiNaC::subs(eval_a, Rs==10.0);
    std::cout << "Ac: " << eval_a.evalf() << std::endl;
    ------------------------------- */








    return 0;
}