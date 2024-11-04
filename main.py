import matplotlib.pyplot as plt
import numpy as np
from numpy import random

import time
import streamlit as st


class Param_Construtivos:
    def __init__(self, Mola, Atrito):
        self.D_ho = 0.0540
        self.D_hp = 0.0635
        self.D_c  = 0.230
        self.ID_s = 0.130
        self.OD_s = 0.160

        self.D_pa = 0.130
        self.R_pa  = 0.065

        self.D_vbhp = 1.2 * self.D_hp
        self.L_vbhp1 = 0.130
        self.L_vbhp2 = 0.010

        self.x_co  = 0.020
        self.x_tot = 0.150000

        self.k_ms = 30000 #kN/m
        self.c_ms = 0.005

        self.PNTA = 5000

        self.L_0 = 0.745200
        self.L_1 = 0.609304
        self.L_2 = 0.458804

        if Mola[0] == True:
            self.k1 = abs(1817000 * Atrito[-1])
            self.k2 = abs(1351000 * Atrito[-1])
        else:
            self.k1 = 1817000
            self.k2 = 1351000


        #====================Simulação de Desgaste do fator de atrito============
        if Atrito[0] == True:
            self.u_rm = abs(0.015 * Atrito[-1])
        else:
            self.u_rm = 0.015

        if Atrito[0] == True:
            self.u_mm = abs(0.120 * Atrito[-1])
        else:
            self.u_mm = 0.120


        #===Dados que o Mashiba passou mais ou menos
        self.D_vpc = 0.105
        self.L_vpc = 0.065

        self.D_vpho  = 0.060
        self.L_vpho1 = 0.035
        self.L_vpho2 = 0.055
        # ===Dados que o Mashiba passou mais ou menos


        self.f_imvpc  = 100#Mashiba não passou...CHUTE!!!
        self.f_imvpho = 100#Mashiba não passou...CHUTE!!!
        self.f_imvbhp = 100#Mashiba não passou...CHUTE!!!


        self.f_imsg = 2 * self.u_mm * self.k_ms * self.c_ms
        self.D_vg = self.ID_s + 1/3*(self.OD_s - self.ID_s)



class Param_Ambientes:
    def __init__(self, Caso):
        self.psi_to_Pa = 6894.76
        self.pa_to_psi = 0.000145038

        self.P_atm = 101325# Pa

        # self.P_tr = 0
        # self.LDA = 2600*0
        # self.AG  = 30*0

        self.D = 0.130
        self.L = 13380
        self.e = 0.00025
        self.v = 2.7*10**-6
        self.API = 25

        # self.P_tr = 10000 * self.psi_to_Pa#psi para Pa
        self.P_sep = 145 * self.psi_to_Pa#psi para Pa

        self.g = 9.806

        self.p_mar  = 1030#998
        self.p_fc   = 1060
        self.p_fcp  = 1083
        self.p_agua = 998

        self.p_fp = (141.5 * self.p_agua)/(self.API + 131.5)

        self.f = 0.017#===CHUTE!!!!!!!!!

        if Caso == 1:
            self.LDA = 0
            self.AG = 0
            self.P_tr = 0
            self.caso = [549.294, 609.273, 609.857, 876.192]
        elif Caso == 2:
            self.LDA = 0
            self.AG = 0
            self.P_tr = 10000 * self.psi_to_Pa  # psi para Pa
            self.caso = [1946, 2006, 1535, 1802]
        elif Caso == 3:
            self.LDA = 2600
            self.AG = 30
            self.P_tr = 0
            self.caso = [89.33, 149.308, 149.893, 416.227]
        elif Caso == 4:
            self.LDA = 2600
            self.AG = 30
            self.P_tr = 10000 * self.psi_to_Pa  # psi para Pa
            self.caso = [1486, 1546, 1076, 1342]
        else:
            self.LDA = 2600
            self.AG = 30
            self.P_tr = 10000 * self.psi_to_Pa  # psi para Pa
            self.caso = [1320, 1380, 1147, 1342]


class Graf_iniciais():
    def __init__(self, Caso):
        self.psi_to_Pa = 6894.76
        self.pa_to_psi = 0.000145038
        self.Caso = Caso

    def graf(self, stroke, Mola, Atrito, ML):
        param_C = Param_Construtivos(Mola, Atrito)
        param_A = Param_Ambientes(Caso=self.Caso)

        F_cm    = []#===Mola
        P_mo    = []#===Pressao Montante
        P_ju    = []#===Pressao Jusante
        P_ju_hx = []#===Pressao Jusante em função da abertura do obturador...
        Kx_x    = []#===resistência ao Fluxo
        vx_x    = []#===Velocidade media do escoamento do fluido

        pos  = []
        por  = []
        hx_x = []

        #===Vetores para guardar os valores
        P_h  = []
        P_sc = []
        P_cfc = []

        F_ph = []
        F_eh = []
        F_sc = []

        f_vpc  = []
        f_vpho = []
        f_vbhp = []
        f_csg  = []
        f_t    = []

        F_at_a = []
        P_at_a = []

        KPs_ar = []

        for ct, i in enumerate(range(0, int(stroke*1000))):
            pos.append(i/1000)

            F_cm.append((-param_C.k1 * (param_C.L_1 - pos[-1])**2) + (param_C.k2 * (param_C.L_1 - pos[-1])))
            P_mo.append(param_A.P_tr)

            #===Porcentagem de passagem no Obturador
            if i <= int(param_C.x_co * 1000):
                P_ju.append(0)
                hx_x.append(0)

                Kx_x.append((1984 * np.exp(-0.735 * (hx_x[-1] ** 0.545))) + 0.1)

                # ===Calculo da velocidade
                p01 = param_A.P_tr - param_A.P_sep - (param_A.p_fp * param_A.g * (param_A.LDA + param_A.AG))
                p02 = param_A.p_fp * param_A.g * (param_A.LDA + param_A.AG)
                p1 = p01 - p02
                p2 = param_A.p_fp*(1 + Kx_x[-1] + (param_A.f*(param_A.L/param_A.D)))
                vx_x.append(np.sqrt(((2 * p1) / p2)))

                # ===Calculo pressão a jusante levando em conta a resistência ao fluxo
                if param_A.P_tr == 0:
                    P_ju_hx.append(0)
                else:
                    P_ju_hx.append(param_A.P_sep + (param_A.p_fp * param_A.g * (param_A.LDA + param_A.AG)))

            else:
                P_ju.append(param_A.P_tr)
                hx_x.append(((pos[-1] - param_C.x_co) / (stroke - param_C.x_co)) * 100)

                Kx_x.append((1984 * np.exp(-0.735 * (hx_x[-1] ** 0.545))) + 0.1)

                # ===Calculo da velocidade
                p01 = param_A.P_tr - param_A.P_sep - (param_A.p_fp * param_A.g * (param_A.LDA + param_A.AG))
                p02 = param_A.p_fp * param_A.g * (param_A.LDA + param_A.AG)
                p1 = p01 + p02
                p2 = param_A.p_fp*(1 + Kx_x[-1] + (param_A.f*(param_A.L/param_A.D)))
                vx_x.append(np.sqrt(((2 * p1)/ p2)))

                #===Calculo pressão a jusante levando em conta a resistência ao fluxo
                if param_A.P_tr == 0:
                    P_ju_hx.append(0)
                else:
                    if self.Caso == 5:
                        P_ju_hx.append(P_mo[-1] - (Kx_x[-1] * param_A.p_fp * ((vx_x[-1]**2)/2)))
                    else:
                        P_ju_hx.append(P_mo[-1])



            #================CALCULOS DO AVANÇO=======================
            #--Iteração Inicial
            if ct == 0:
                P_at_aux = (param_A.p_mar * param_A.g * param_A.LDA) / (np.pi/4 * (param_C.D_c**2 - param_C.D_ho**2))
            else:
                P_at_aux = P_at_a[-1]

            #---Pressão Hidrostática devido a LDA
            P_h.append((param_A.p_mar * param_A.g * param_A.LDA)) #--> Pa
            P_cfc.append(param_A.p_fc * param_A.g * (param_A.LDA + param_A.AG)) #--> Pa

            #---Pressão->Força Hidrostática na haste que fica exposta ao mar (Overrride)
            F_ph.append(P_h[-1] * ((np.pi/4*(param_C.D_ho**2))))

            #---Pressão->Força da Cavidade na Ponta da Haste
            F_eh.append(param_A.P_tr * ((np.pi/4*(param_C.D_hp**2))))

            #---Pressão->Força do Sistema de Compensão na Coroa Circular "Interna"
            F_sc.append(P_h[-1] * ((np.pi/4*(param_C.D_c**2 - param_C.D_hp**2))))

            # ...Atrito vedação pistão-cilindro
            f_vpc.append((param_C.f_imvpc + param_C.u_rm * (P_at_aux - P_h[-1]) * np.pi * param_C.D_vpc * param_C.L_vpc)*0)

            #...Atrito pistão-haste-override
            f_vpho.append((param_C.f_imvpho + param_C.u_rm * np.pi * param_C.D_vpho *
                      ((P_at_aux - param_A.P_atm)*param_C.L_vpho1 + (P_h[-1] - param_A.P_atm)*param_C.L_vpho2))*0)

            #...Atrito vedação bonnet-haste principal
            f_vbhp.append((param_C.f_imvbhp + param_C.u_rm * np.pi * param_C.D_vbhp *
                       (P_at_aux - param_A.P_atm) * param_C.L_vbhp1 +
                       (P_cfc[-1] - param_A.P_atm) * param_C.L_vbhp2)*0)

            #...Atrito contato sede-gaveta
            f_csg.append((param_C.f_imsg + (np.pi/4 * param_C.u_mm * (P_mo[-1] - P_ju_hx[-1]) * param_C.D_vg**2))*1)

            #===Atrito TOTAL
            f_t.append(f_csg[-1] + f_vpc[-1] + f_vpho[-1] + f_vbhp[-1])

            #===Força de Atuação TOTAL
            F_at_a.append(-F_ph[-1] + F_eh[-1] + F_sc[-1] + f_t[-1] + F_cm[-1])

            #===Pressão de Atuação - P_at_av
            if ct == 0:
                P_at_a.append(P_h[0])
            elif ct == int(stroke*1000)-1:
                P_at_a.append((param_C.PNTA*param_A.psi_to_Pa) + P_cfc[0])
            else:
                P_at_a.append((F_at_a[-1] / (np.pi/4 * (param_C.D_c**2 - param_C.D_ho**2))))


            # if hx_x[-1] > 29:
            #     if hx_x[-1] < 31:
            #         print('Porcentagem {} com perda de carga de {}'.format(hx_x[-1], Kx_x[-1]))


        fig1 = plt.figure()

        #===MOLA==========================
        ax = fig1.add_subplot(2, 3, 1)
        ax.plot(pos, F_cm, c='r')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Força Mola [N]")
        ax.set_title("Força Mola Belleville")
        ax.grid()
        ax1 = ax.twinx()
        ax1.plot(pos, hx_x, c='b', ls='--', lw='0.45')
        ax1.set_ylabel("h(x)")

        #===JUSANTE========================
        ax = fig1.add_subplot(2, 3, 2)
        ax.plot(pos, P_mo, c='green', ls='--', label='Montante')
        ax.plot(pos, P_ju, c='b'    , ls='-', label='Jusante')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Pressão Mon x Jus [Pa]")
        ax.set_title("Montante e Jusante")
        ax.legend()
        ax1 = ax.twinx()
        ax1.plot(pos, hx_x, c='b', ls='--', lw='0.45')
        ax1.set_ylabel("h(x)")
        ax.grid()


        #===RESISTENCIA AO FLUXO============
        ax = fig1.add_subplot(2, 3, 3)
        ax.plot(pos, Kx_x, c='r', ls='-')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Resistência ao Fluxo [-]")
        ax.set_title("Resistência ao Escoamento")
        ax.grid()
        ax1 = ax.twinx()
        ax1.plot(pos, hx_x, c='b', ls='--', lw='0.45')
        ax1.set_ylabel("h(x)")

        #===VELOCIDADE DO ESCOAMENTO============
        ax = fig1.add_subplot(2, 3, 4)
        ax.plot(pos, vx_x, c='b'    , ls='--', label='Velocidade')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Velocidade [m/s]")
        ax.set_title("Velocidade")
        ax.legend()
        ax.grid()
        ax1 = ax.twinx()
        ax1.plot(pos, hx_x, c='b', ls='--', lw='0.45')
        ax1.set_ylabel("h(x)")


        #===PRESSÃO NA JUSANTE CORRIGIDO============
        ax = fig1.add_subplot(2, 3, 5)
        ax.plot(pos, P_ju_hx, c='r', ls='-', label='Jusante Corrigida')

        ax.plot(pos, P_mo, c='green', ls='--', label='Montante')
        ax.plot(pos, P_ju, c='b'    , ls='--', label='Jusante')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Resistência ao Fluxo [-]")
        ax.set_title("Jusante Corrigida")
        ax.legend()
        ax.grid()
        ax1 = ax.twinx()
        ax1.plot(pos, hx_x, c='b', ls='--', lw='0.45')
        ax1.set_ylabel("h(x)")

        # fig1.tight_layout()

        # ===Assinatura==========================
        fig2 = plt.figure()
        ax = fig2.add_subplot(2, 2, 1)
        ax.plot(pos, P_at_a, c='b', lw=0.75, label='pressão')
        ax.plot(pos, P_cfc, c='g',  ls='-', lw=0.99, label='Sistema Compensação')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Pressão de Atuação [Pa]")
        ax.legend()
        ax.grid()

        # ax = fig2.add_subplot(2, 2, 2)
        # ax.plot(pos, [n * 0.000145038 for n in P_at_a], c='r', lw=0.75, label='pressão')
        # ax.plot(pos, [n * 0.000145038 for n in P_cfc], c='g', ls='-', lw=0.99, label='Sistema Compensação')
        # ax.set_xlabel("Posição Haste (Stroke) [m]")
        # ax.set_ylabel("Pressão de Atuação [psi]")
        # ax.set_title("Assinatura Mashiba")
        # ax.legend()
        # ax.grid()

        a = list(np.array(P_at_a) - np.array(P_cfc))
        pre_psi = [n * 0.000145038 for n in a]
        ax = fig2.add_subplot(2, 2, 2)
        ax.plot(pos, [n * 0.000145038 for n in a], c='r', lw=0.75, label='Pressão (- P_cfc)')
        ax.plot(pos, [n * 0.000145038 for n in P_cfc], c='g', ls='-', lw=0.99, label='Sistema Compensação')
        ax.set_xlabel("Posição Haste (Stroke) [m]")
        ax.set_ylabel("Pressão de Atuação [psi]")
        ax.legend()
        ax.grid()


        #====Key Points=========
        b = [n * 0.000145038 for n in a]
        ax = fig2.add_subplot(2, 2, 3)
        posicoes = [1, (param_C.x_co*1000-1), (param_C.x_co*1000+1), (param_C.x_tot*1000-2)]

        dif = []
        KPs = ['A2',  'A3', 'A4', 'A5']
        kp_aux = []
        for ct, (kp, i, j) in enumerate(zip(KPs, param_A.caso, posicoes)):
            ax.scatter(kp, i, color='green')
            ax.scatter(kp, b[int(j)], color='red')
            dif.append(i - b[int(j)])
            kp_aux.append(b[int(j)])
        KPs_ar.append(kp_aux)
        ax.set_xlabel("KP [-]")
        ax.set_ylabel("Pressão de Atuação [psi]")
        ax.legend()
        ax.grid()


        ax = fig2.add_subplot(2, 2, 4)
        ax.bar(KPs, dif)
        ax.set_xlabel("KP [-]")
        ax.set_ylabel("Dif entre KPs [psi]")
        ax.legend()
        ax.grid()

        fig2.suptitle("Caso: " + str(self.Caso))

        if ML == 0:
            nada = 0
        else:
            plt.close(fig1)
            plt.close(fig2)
        return KPs_ar, pre_psi, dif, F_cm, Kx_x, vx_x



def chama_funcoes(atrito):
    comp_1 = []
    comp_2 = []
    Estudar = ['Mola', 'Atrito']#===Gera os parâmetros que serão estudados. ex: Mola e Atrito
    Estudar = ['Nada']#===Gera os parâmetros que serão estudados. ex: Mola e Atrito

    for cond in Estudar:
        #===Cria vetores para armazenar as variaveis
        A2 = []
        A3 = []
        A4 = []
        A5 = []
        Ass = []

        #===Cria os plots e subplots
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)

        #===Cria
        fig21 = plt.figure()
        ax_21 = fig21.add_subplot(2, 1, 1)

        a = Graf_iniciais(Caso=4)#====Define o Caso estudado pelo MASHIBA

        ML = False

        for i in range(1):#===Gera os cenarios de estudo========================================
            # val_Atrito = random.normal(loc=0, scale=1)
            val_Atrito = atrito

            Bol_Mola = True
            Bol_Atri = True

            if cond == 'Mola':
                Bol_Mola = True
            if cond == 'Atrito':
                Bol_Atri = True


            KPs_ar, P_at_a, dif_KPs, F_cm, Kx_x, vx_x = a.graf(stroke=0.150, Mola=[Bol_Mola, val_Atrito], Atrito=[Bol_Atri, abs(val_Atrito)], ML=i)
            ax_21.plot(P_at_a)

            Ass.append(P_at_a)
            A2.append(KPs_ar[0][0])
            A3.append(KPs_ar[0][1])
            A4.append(KPs_ar[0][2])
            A5.append(KPs_ar[0][3])

            ax1.scatter(A2[-1], A3[-1], edgecolors='k', s=100)
            ax2.scatter(A2[-1], A4[-1], edgecolors='k', s=100)
            ax3.scatter(A2[-1], A5[-1], edgecolors='k', s=100)
            ax4.scatter(A3[-1], A4[-1], edgecolors='k', s=100)

            ax1.set_xlabel("A2")
            ax1.set_ylabel("A3")

            ax2.set_xlabel("A2")
            ax2.set_ylabel("A4")

            ax3.set_xlabel("A2")
            ax3.set_ylabel("A5")

            ax4.set_xlabel("A3")
            ax4.set_ylabel("A4")

            fig.suptitle(cond)

        # if ML == True:
    #         from sklearn.decomposition import PCA
    #         from sklearn.preprocessing import MinMaxScaler, StandardScaler, Normalizer
    #
    #         ass = np.array(Ass)
    #         scale = MinMaxScaler()
    #         # scale = StandardScaler()
    #         # scale = Normalizer()
    #         scaledX = scale.fit_transform(ass.T)
    #
    #         pca = PCA(n_components=3)
    #         X_train = pca.fit_transform(np.array(scaledX.T))
    #
    #         fig = plt.figure()
    #         ax1 = fig.add_subplot(2, 2, 1)
    #         ax1.scatter(X_train[:, 0], X_train[:, 1], edgecolors='k')
    #         ax1.set_xlabel("Comp #1")
    #         ax1.set_ylabel("Comp #2")
    #         fig.suptitle(cond)
    #
    #         comp_1.append(X_train[:, 0])
    #         comp_2.append(X_train[:, 1])
    #
    # if ML == True:
    #     fig = plt.figure()
    #     ax = fig.add_subplot(1, 1, 1)
    #     for i in range(len(comp_1)):
    #         ax.scatter(comp_1[i], comp_2[i], edgecolors='k')

    plt.show()

    return KPs_ar, P_at_a, dif_KPs, F_cm, Kx_x, vx_x

def main():
    st.title('Assinatura Avanço Mashiba')
    st.text('Simulação da assinatura de pressão durante o avanço do pistão')

    st.header("Valor de Modificação")
    st.text("Atrito e Mola Simultâneos...")
    slider_value1 = st.slider("Escolha um valor", 0, 100, 1) / 100
    st.write("Valor escolhido: ", slider_value1)

    KPs_ar, P_at_a, dif_KPs, F_cm, Kx_x, vx_x = chama_funcoes(slider_value1)

    st.header("Gráfico Assinatura LAMEF")
    st.line_chart(P_at_a, x_label="Stroke (mm)", y_label="Pressão (psi)")

    st.header("Força Mola LAMEF")
    st.line_chart(F_cm, x_label="Stroke (mm)", y_label="Força Mola (N)")

    st.header("Resistência ao Escoamento LAMEF")
    st.line_chart(Kx_x, x_label="Stroke (mm)", y_label="Resis. Fluxo (-)")

    st.header("Velocidade LAMEF")
    st.line_chart(vx_x, x_label="Stroke (mm)", y_label="Velocidade (m/s)")

    st.header("Diferença entre modelo LAMEF e MASHIBA")
    data = {"x": ['A2', 'A3', 'A4', 'A5'], "y": dif_KPs}
    st.bar_chart(data, x='x', y='y', x_label="Key Points (KPs)", y_label="Pressão (psi)")

    return KPs_ar, P_at_a, dif_KPs

main()