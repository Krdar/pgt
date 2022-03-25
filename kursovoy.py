from iapws import IAPWS97 as WSP
from iapws import IAPWS97 as IAPWS97 
import math as M
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import pandas as pd

st.title('Курсовая работа часть 2')
st.write('Расчет регулирующей ступени паровой турбины Т-180/215-12,8-2 ЛМЗ')
st.subheader('Крыницкая Д. ФПэ-01-19')


def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color = 'с'):
    if not isinstance(wsp_point,list):
        iso_bar_0_s = np.arange(wsp_point.s+min_s,wsp_point.s+max_s,step_s).tolist()
        iso_bar_0_h = [WSP(P = wsp_point.P, s = i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s+min_s,wsp_point[1].s+max_s,step_s).tolist()
        iso_bar_0_h = [WSP(P = wsp_point[1].P, s = i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s,iso_bar_0_h,color)



st.write("Давление пара перед ступенью:")
st.write("p_0 = 13.5e6")
st.write("Температура пара перед ступенью:")
st.write("T_0 = 552 + 273.15")
st.write("Расход пара через ступень:")
st.write("G_0 = 177.413")
st.write("Располагаемый теплоперепад:")
st.write("H0= 90-110")
st.write("Скорость пара на входе в сопловую решетку:")
st.write("c0 = 0")
st.write("Степень реактивности:")
st.write("ρ=0.05")
st.write("Частота вращения ротора турбины:")
st.write("n = 50 ")
st.write("Хорда сопловой решетки:")
st.write("b_1 = 0.06")
st.write(" Средний диаметр регулирующей ступени: ")
st.write("d = 1.1")
st.write("Хорда рабочей решетки:")
st.write("b_2 = 0.03")
st.write("delta = 0.003")
st.write("Высота сопловых лопаток:")
st.write("l_1 = 0.015")
st.write("Угол выхода из сопловой решетки:")
st.write("alpha_1e = 14")
st.write("Коэффициент использования выходной скорости: ")
st.write("kappa_vs=0")

p_0 = 13.5
T_0 = 552 + 273.15
G_0 = 177.413
H_0= 90
c0 = 0
rho=0.05
n = 50 
b_1 = 0.06
d = 1.1
b_2 = 0.03

l_1 = 0.015
alpha_1 = 14
Delta = 0.003 
kappa_vs=0

st.write(" ")
st.write("**Решение:** ")

def callculate_optimum(d, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs):
    u = M.pi*d*n
    point_0 = WSP(P = p_0, T = T_0)
    H_0s = H_0*(1-rho)
    H_0r = H_0*rho
    h_1t = point_0.h - H_0s
    point_1t = WSP(h = h_1t, s = point_0.s)
    c_1t = (2000*H_0s)**0.5
    M_1t = c_1t/point_1t.w
    mu_1 = 0.982 - 0.005*(b_1/l_1)
    F_1 = G_0*point_1t.v/mu_1/c_1t
    el_1 = F_1/M.pi/d/M.sin(M.radians(alpha_1))
    e_opt=5*el_1**0.5
    if e_opt > 0.85:
        e_opt = 0.85
    l_1 = el_1/e_opt
    fi_1 = 0.98 - 0.008 * (b_1 / l_1)
    c_1 = c_1t * fi_1
    alpha_1 = M.degrees(M.asin(mu_1 / fi_1 * M.sin(M.radians(alpha_1))))
    w_1 = (c_1 ** 2 + u ** 2 - 2 * c_1 * u * M.cos(M.radians(alpha_1))) ** 0.5
    betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1)) / (M.cos(M.radians(alpha_1)) - u / c_1)))
    Delta_Hs = c_1t ** 2 / 2 * (1 - fi_1 ** 2)
    h_1 = h_1t + Delta_Hs * 1e-3
    point_1 = WSP(P=point_1t.P, h=h_1)
    h_2t = h_1 - H_0r
    point_2t = WSP(h=h_2t, s=point_1.s)
    w_2t = (2 * H_0r * 1e3 + w_1 ** 2) ** 0.5
    l_2 = l_1 + Delta
    mu_2 = 0.965 - 0.01 * (b_2 / l_2)
    M_2t = w_2t / point_2t.w
    F_2 = G_0 * point_2t.v / mu_2 / w_2t
    betta_2 = M.degrees(M.asin(F_2 / (e_opt * M.pi * d * l_2)))
    point_1w = WSP(h=point_1.h + w_1 ** 2 / 2 * 1e-3, s=point_1.s)
    psi = 0.96 - 0.014 * (b_2 / l_2)
    w_2 = psi * w_2t
    c_2 = (w_2 ** 2 + u ** 2 - 2 * u * w_2 * M.cos(M.radians(betta_2))) ** 0.5
    alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2)) / (M.cos(M.radians(betta_2)) - u / w_2)))
    if alpha_2 < 0:
        alpha_2 = 180 + alpha_2
    Delta_Hr = w_2t ** 2 / 2 * (1 - psi ** 2)
    h_2 = h_2t + Delta_Hr * 1e-3
    point_2 = WSP(P=point_2t.P, h=h_2)
    Delta_Hvs = c_2 ** 2 / 2
    E_0 = H_0 - kappa_vs * Delta_Hvs
    etta_ol1 = (E_0 * 1e3 - Delta_Hs - Delta_Hr - (1 - kappa_vs) * Delta_Hvs) / (E_0 * 1e3)
    etta_ol2 = (u * (c_1 * M.cos(M.radians(alpha_1)) + c_2 * M.cos(M.radians(alpha_2)))) / (E_0 * 1e3)
    return etta_ol2, alpha_2

st.write(" ")
d = [i * 1e-2 for i in list(range(90, 111, 1))]

eta = []
ucf = []
al=[]
for i in d:
    ucf_1 = M.pi * i * n / (2000 * H_0) ** 0.5
    ucf.append(ucf_1)
    eta_ol, alpha = callculate_optimum(i, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs)
    eta.append(eta_ol)
    al.append(alpha)
    
eta_=np.around(eta, decimals=3)
ucf_=np.around(ucf, decimals=3)
al_=np.around(al, decimals=2)

ucf_eta = plt.figure()

plt.plot(ucf, eta)
plt.title("Зависимость ηo.л от u/cf")
plt.xlabel("u/cf")
plt.ylabel("ηo.л")
plt.grid()

st.pyplot(ucf_eta)
w=[i for i in list(range(90, 111, 1))]
import pandas as pd
itog = pd.DataFrame ({
                      'ηо.л':eta_,
                      '(u/сф)':ucf_,
                      'α2':al_,
                       'H_0, кДж/кг': w,
                     })
itog

st.write(" ")
d = 1.1
u = M.pi*d*n
print(f'u = {u:.2f} м/с')
point_0 = WSP(P = p_0, T = T_0)
print(f'h_0 = {point_0.h:.2f} кДж/кг')
print(f's_0 = {point_0.s:.4f} кДж/(кг*К)')
H_0s = H_0*(1-rho)
H_0r = H_0*rho
h_1t = point_0.h - H_0s
print(f'h_1т = {h_1t:.2f} кДж/кг')
point_1t = WSP(h = h_1t, s = point_0.s)
c_1t = (2000*H_0s)**0.5
print(f'c_1т = {c_1t:.2f} м/с')
M_1t = c_1t/point_1t.w
print(f'M_1т = {M_1t:.2f}')
mu_1 = 0.982 - 0.005*(b_1/l_1)
F_1 = G_0*point_1t.v/mu_1/c_1t
print(f'F_1 = {F_1:.4f} м^2')
el_1 = F_1/M.pi/d/M.sin(M.radians(alpha_1))
print(f'el_1 = {el_1:.4f} м')
e_opt=6*el_1**0.5
if e_opt > 0.85:
    e_opt = 0.85
l_1 = el_1/e_opt

hs = plt.figure()

def plot_hs_nozzle_t(x_lim, y_lim):
    plt.plot([point_0.s, point_1t.s],[point_0.h, point_1t.h],'ro-')
    iso_bar(point_0,-0.02,0.02,0.001,'c')
    iso_bar(point_1t,-0.02,0.02,0.001,'m')
    plt.xlim(x_lim)
    plt.ylim(y_lim)

plot_hs_nozzle_t([6.5, 6.7],[3300,3500])
plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")
plt.grid()

#st.pyplot(hs)

st.write(" ")
print(f'l_1 = {l_1:.4f} м')
if alpha_1 <= 10:
    NozzleBlade = 'C-90-09A'
    t1_ = 0.78
    b1_mod = 6.06
    f1_mod = 3.45
    W1_mod = 0.471
    alpha_inst1 = alpha_1-12.5*(t1_-0.75)+20.2
elif  10 < alpha_1 <= 13:
    NozzleBlade = 'C-90-12A'
    t1_ = 0.78
    b1_mod = 5.25
    f1_mod = 4.09
    W1_mod = 0.575
    alpha_inst1 = alpha_1-10*(t1_-0.75)+21.2
elif  13 < alpha_1 <= 16:
    NozzleBlade = 'C-90-15A'
    t1_ = 0.78
    b1_mod = 5.15
    f1_mod = 3.3
    W1_mod = 0.45
    alpha_inst1 = alpha_1-16*(t1_-0.75)+23.1
else:
    NozzleBlade = 'C-90-18A'
    t1_ = 0.75
    b1_mod = 4.71
    f1_mod = 2.72
    W1_mod = 0.333
    alpha_inst1 = alpha_1-17.7*(t1_-0.75)+24.2
print('*Тип профиля:*',  NozzleBlade)
print(f'*Оптимальный относительный шаг t1_* = {t1_}')
z1 = (M.pi*d)/(b_1*t1_)
z1 = int(z1)
if z1 % 2 == 0:
    print(f'z1 = {z1}')
else:
    z1 = z1+1
    print(f'z1 = {z1}')
t1_ = (M.pi*d)/(b_1*z1)
Ksi_1_ = (0.021042*b_1/l_1 + 0.023345)*100
k_11 = 7.18977510*M_1t**5 - 26.94497258*M_1t**4 + 39.35681781*M_1t**3 - 26.09044664*M_1t**2 + 6.75424811*M_1t + 0.69896998
k_12 = 0.00014166*90**2 - 0.03022881*90 + 2.61549380
k_13 = 13.25474043*t1_**2 - 20.75439502*t1_ + 9.12762245
Ksi_1 = Ksi_1_*k_11*k_12*k_13

fi_1 = M.sqrt(1-Ksi_1/100)
print(f'mu_1 = {mu_1}')
print(f'fi_1 = {fi_1}')

st.write(" ")
alpha_1 = 12
c_1 = c_1t*fi_1
print(f'c_1 = {c_1:.2f} м/с')
alpha_1 = M.degrees(M.asin(mu_1/fi_1*M.sin(M.radians(alpha_1))))
print(f'alpha_1 = {alpha_1:.2f} град.')
w_1 = (c_1**2+u**2-2*c_1*u*M.cos(M.radians(alpha_1)))**0.5
print(f'w_1 = {w_1:.4f}')

cwu = plt.figure()

c_1u = c_1*M.cos(M.radians(alpha_1))
c_1a = c_1*M.sin(M.radians(alpha_1))
w_1u = c_1u - u
print(f'{c_1u:.4f}', f' ;  {w_1u:.4f}')
w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]
ax = plt.axes()
ax.arrow(*c_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
plt.text(-2*c_1u/3, -3*c_1a/4, '$c_1$', fontsize=20)
plt.text(-2*w_1u/3, -3*c_1a/4, '$w_1$', fontsize=20)
plt.title("Входной треугольник скоростей")
plt.grid()

#st.pyplot(cwu)


st.write(" ")
betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1))/(M.cos(M.radians(alpha_1))-u/c_1)))
print(f'betta_1 = {betta_1:.2f}')
Delta_Hs = c_1t**2/2*(1-fi_1**2)
print(f'Delta_Hs = {Delta_Hs:.2f} Дж/кг')

st.write(" ")
h_1 = h_1t + Delta_Hs*1e-3
print(f'h_1 = {h_1:.2f} кДж/кг')
point_1 = WSP(P = point_1t.P, h = h_1)
h_2t = h_1 - H_0r
print(f'h_2t = {h_2t:.2f} кДж/кг')
point_2t = WSP(h = h_2t, s = point_1.s)
w_2t = (2*H_0r*1e3+w_1**2)**0.5
print(f'w_2t = {w_2t:.2f} м/с')
l_2 = l_1 + Delta
mu_2 = 0.965-0.01*(b_2/l_2)
print(f'mu_2 = {mu_2:.2f}')
M_2t = w_2t/point_2t.w
print(f'M_2t = {M_2t:.2f}')
F_2 = G_0*point_2t.v/mu_2/w_2t
print(f'F_2 = {F_2:.2f}')
betta_2 = M.degrees(M.asin(F_2/(e_opt*M.pi*d*l_2)))
print(f'betta_2 = {betta_2:.2f}')
point_1w = WSP(h = point_1.h+w_1**2/2*1e-3, s = point_1.s)
fig4 = plt.figure(figsize=(10, 10))
psi = 0.91

w_2 = psi*w_2t
print(f'w_2 = {w_2:.2f} м/с')
c_2 = (w_2**2+u**2-2*u*w_2*M.cos(M.radians(betta_2)))**0.5
print(f'c_2 = {c_2:.2f} м/с')
alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2))/(M.cos(M.radians(betta_2))-u/w_2)))
print(f'alpha_2 = {alpha_2:.2f}')
Delta_Hr = w_2t**2/2*(1-psi**2)
print(f'Delta_Hr = {Delta_Hr:.2f} Дж/кг')
h_2 = h_2t+Delta_Hr*1e-3
point_2 = WSP(P = point_2t.P, h = h_2)
Delta_Hvs = c_2**2/2
print(f'Delta_Hvs = {Delta_Hvs:.2f} Дж/кг')
E_0 = H_0 - kappa_vs*Delta_Hvs
etta_ol1 = (E_0*1e3 - Delta_Hs-Delta_Hr-(1-kappa_vs)*Delta_Hvs)/(E_0*1e3)
print(f'1. etta_ol = {etta_ol1}')
etta_ol2 = (u*(c_1*M.cos(M.radians(alpha_1))+c_2*M.cos(M.radians(alpha_2))))/(E_0*1e3)
print(f'2. etta_ol = {etta_ol2}')


h_3 = h_2 + Delta_Hvs * 1e-3
point_3 = IAPWS97(P=point_2t.P, h=h_3)

point_2_ = IAPWS97(P=point_2t.P, h=point_0.h-H_0)

hsstage = plt.figure()

def plot_hs_stage_t(x_lim,y_lim):
    plot_hs_nozzle_t(x_lim,y_lim)
    plt.plot([point_0.s,point_1.s],[point_0.h,point_1.h],'bo-')
    plt.plot([point_1.s,point_2t.s],[point_1.h,point_2t.h], 'ro-')
    plt.plot([point_1.s,point_1.s],[point_1w.h, point_1.h],'ro-')
    plt.plot([point_1.s, point_2.s], [point_1.h, point_2.h], 'bo-')
    plt.plot([point_2.s, point_3.s], [point_2.h, point_3.h], 'bo-')
    iso_bar(point_2t,-0.02,0.02,0.001,'y')
    iso_bar(point_1w,-0.005,0.005,0.001,'c')
plot_hs_stage_t([6.55,6.65],[3300,3500])
plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")
plt.grid()

st.pyplot(hsstage)

st.write(" ")
if betta_2 <= 15:
    RotorBlade = 'P-23-14A'
    t2_ = 0.63
    b2_mod = 2.59
    f2_mod = 2.44
    W2_mod = 0.39
    beta_inst2 = betta_2-12.5*(t2_-0.75)+20.2

elif  15 < betta_2 <= 19:
    RotorBlade = 'P-26-17A'
    t2_ = 0.65
    b2_mod = 2.57
    f2_mod = 2.07
    W2_mod = 0.225
    beta_inst2 = betta_2-19.3*(t2_-0.6)+60


elif  19 < betta_2 <= 23:
    RotorBlade = 'P-30-21A'
    t2_ = 0.63
    b2_mod = 2.56
    f2_mod = 1.85
    W2_mod = 0.234
    beta_inst2 = betta_2-12.8*(t2_-0.65)+58


elif 23 < betta_2 <= 27:
    RotorBlade = 'P-35-25A'
    t2_ = 0.6
    b2_mod = 2.54
    f2_mod = 1.62
    W2_mod = 0.168
    beta_inst2 = betta_2-16.6*(t2_-0.65)+54.3

elif 27 < betta_2 <= 31:
    RotorBlade = 'P-46-29A'
    t2_ = 0.51
    b2_mod = 2.56
    f2_mod = 1.22
    W2_mod = 0.112
    beta_inst2 = betta_2-50.5*(t2_-0.6)+47.1


else:
    RotorBlade = 'P-50-33A'
    t2_ = 0.49
    b2_mod = 2.56
    f2_mod = 1.02
    W2_mod = 0.079
    beta_inst2 = betta_2-20.8*(t2_-0.6)+43.7

print('*Тип профиля:*',  RotorBlade)
print(f'*Оптимальный относительный шаг t2_* = {t2_}')

z2 = int((M.pi*d)/(b_2*t2_))
print(f'z2 = {z2}')
t2_ = (M.pi*d)/(b_2*z2)
Ksi_2_ = 4.364*b_2/l_2 + 4.22
k_21 = -13.79438991*M_2t**4 + 36.69102267*M_2t**3 - 32.78234341*M_2t**2 + 10.61998662*M_2t + 0.28528786
k_22 = 0.00331504*betta_1**2 - 0.21323910*betta_1 + 4.43127194
k_23 = 60.72813684*t2_**2 - 76.38053189*t2_ + 24.97876023
Ksi_2 = Ksi_2_*k_21*k_22*k_23

psi = M.sqrt(1-Ksi_2/100)
print(f'psi = {psi:.2f}')



st.write(" ")
cw = plt.figure()

c_1u = c_1*M.cos(M.radians(alpha_1))
c_1a = c_1*M.sin(M.radians(alpha_1))
w_1u = c_1u - u
w_2a = w_2*M.sin(M.radians(betta_2))
w_2u = w_2*M.cos(M.radians(betta_2))
c_2u=w_2u-u
print(f'c_1u = {c_1u}')
print(f'w_1u = {w_1u}')
w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]

w_2_tr = [0, 0, w_2u, -w_2a]
c_2_tr = [0, 0, c_2u, -w_2a]
u_2_tr = [w_2u,-w_2a, -u, 0]
ax = plt.axes()
ax.arrow(*c_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
ax.arrow(*c_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
plt.text(-2*c_1u/3, -3*c_1a/4, '$c_1$', fontsize=20)
plt.text(-2*w_1u/3, -3*c_1a/4, '$w_1$', fontsize=20)
plt.text(2*c_2u/3, -3*w_2a/4, '$c_2$', fontsize=20)
plt.text(2*w_2u/3, -3*w_2a/4, '$w_2$', fontsize=20)
plt.title("Входной и выходной треугольники скоростей")
plt.grid()

st.pyplot(cw)

st.write(" ")
delta_a = 0.0025
z_per_up = 2
mu_a = 0.5
mu_r = 0.75
d_per = d + l_1
delta_r = d_per*0.001
delta_ekv = 1/M.sqrt(1/(mu_a*delta_a)**2+z_per_up/(mu_r*delta_r)**2)
print("""Эквивалентный зазор в уплотнении по бандажу (периферийном) 
            delta_ekv = %.3f мм""" % (delta_ekv*1000))
xi_u_b=M.pi*d_per*delta_ekv*etta_ol1/F_1*M.sqrt(rho+1.8*l_2/d)
print("""Относительные потери от утечек через бандажные уплотнения 
                xi_u_b = %.3f""" % xi_u_b)
Delta_Hub = xi_u_b*E_0
print("""Абсолютные потери от утечек через периферийное уплотнение ступени  
            Delta_Hub = %.3f кДж/кг""" % Delta_Hub)

k_tr=0.0007
Kappa_VS = 0
u = M.pi*d*n
c_f = M.sqrt(2000*H_0)
ucf = u/c_f
xi_tr = k_tr*d**2/F_1*ucf**3
print("""Определяем u/c_ф для ступени  
            U/c_ф = %.3f""" % ucf)
print("""Относительные потери от трения диска  
            xi_tr =  = %.5f""" % xi_tr)
Delta_Htr = xi_tr*E_0
print("""Абсолютные потери от трения диска  
            Delta_Htr = %.3f кДж/кг""" % Delta_Htr)

k_v = 0.065
m = 1
xi_v = k_v/M.sin(M.radians(alpha_1))*(1-e_opt)/e_opt*ucf**3*m
print("""Относительные вентиляционные потери 
            xi_v = %.5f """ % xi_v)
i_p = 4
B_2 = b_2*M.sin(M.radians(beta_inst2))
xi_segm = 0.25*B_2*l_2/F_1*ucf*etta_ol1*i_p
print("""Относительные сегментные потери
            xi_segm = %.5f""" % xi_segm)
xi_parc = xi_v+xi_segm
Delta_H_parc = E_0*xi_parc

H_i = E_0 - Delta_Hr*1e-3 - Delta_Hs*1e-3 - (1-Kappa_VS)*Delta_Hvs*1e-3 - Delta_Hub - Delta_Htr - Delta_H_parc
print("""Использованный теплоперепад ступени  
           H_i = %.3f кДж/кг$""" % H_i)

eta_oi = H_i/E_0
st.write("""Внутренний относительный КПД ступени  
        eta_oi  = %.3f """ % eta_oi)
N_i = G_0*H_i
st.write("""Внутреняя мощность ступени  
            N_i = %.2f кВт""" % N_i)
