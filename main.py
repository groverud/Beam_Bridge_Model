'''

||||||||||||||||||||                                             | n_flanges_top
||||||||||||||||||||                                             | h_flange
    ||||    ||||
    ||||    ||||        ------------------- centroidal axis                         | n_webs
    ||||    ||||                                                                    | h_web
    ||||    ||||                                                                    | b_webs
||||||||||||||||||||    __________________ datum                 | n_flanges_bot


'''



import math as m


class Beam:
    def __init__(self, n_webs = 0, n_flanges_top = 0, n_flanges_bot = 0, b_web = 0, b_flange = 0, h_web = 0, h_flange = 0):
        self.n_webs = n_webs
        self.n_flanges_top = n_flanges_top
        self.n_flanges_bot = n_flanges_bot
        self.b_web = b_web
        self.b_flange = b_flange
        self.h_web = h_web
        self.h_flange = h_flange

    def calc_Centroid(self):

        bot_A = self.n_flanges_bot * self.h_flange * self.b_flange
        bot_y = (self.n_flanges_bot* self.h_flange)/2

        mid_A = self.n_webs * self.h_web * self.b_web
        mid_y = (self.h_web/2) + (self.n_flanges_bot * self.h_flange)

        top_A = self.n_flanges_top * self.h_flange * self.b_flange 
        top_y = (self.n_flanges_top * self.h_flange)/2 + (self.h_web) + (self.n_flanges_bot * self.h_flange)

        Y = (top_A * top_y + mid_A * mid_y + bot_A * bot_y)/ (top_A + mid_A + bot_A)

        return Y

    def calc_I(self):

        Y = self.calc_Centroid()

        bot_A = self.n_flanges_bot * self.h_flange * self.b_flange
        bot_y_b = abs(Y - ( (self.n_flanges_bot * self.h_flange)/2 ))
        #After || axis
        bot_I = (self.b_flange * (self.n_flanges_bot * self.h_flange)**3)/12 + (bot_A * (bot_y_b**2))           

        mid_A = self.n_webs * self.h_web * self.b_web
        mid_y_b = abs( Y - ((self.h_web/2) + (self.n_flanges_bot * self.h_flange)) ) 
        #After || axis
        mid_I = self.n_webs*(self.b_web * (self.h_web)**3)/12 + (mid_A * (mid_y_b**2))       

        top_A = self.n_flanges_top * self.h_flange * self.b_flange 
        top_y_b = abs(Y - (self.n_flanges_top * self.h_flange)/2 + (self.h_web) + (self.n_flanges_bot * self.h_flange))
        #After || axis
        top_I = (self.b_flange * (self.n_flanges_top * self.h_flange)**3)/12 + (top_A * (top_y_b**2))    

        I = bot_I + mid_I + top_I  

        return I

    def calc_Q(self):
        Y = self.calc_Centroid()

        Q_glue = (self.n_flanges_top* self.h_flange * self.b_flange)  *  abs(Y - (self.n_flanges_top * self.h_flange)/2 + (self.h_web) + (self.n_flanges_bot * self.h_flange))
        Q_section = (self.n_flanges_top* self.h_flange * self.b_flange + (self.n_webs * self.b_web * self.h_web/2) )  *  abs( Y - ((self.h_web/2) + (self.n_flanges_bot * self.h_flange)) )
        
        return (Q_glue, Q_section)

    
    #Failure P Checks: 

    def calc_P_Flex(self):
        '''
        check for (t and c) at both (190 and 166) mark i.e. 4 values
        '''

        I = self.calc_I()
        Y = self.calc_Centroid()

        Pc_166 = (6 * I)/ (166 * ((self.n_flanges_bot * self.h_flange + self.h_web + self.n_flanges_top * self.h_flange)-Y))
        Pt_166 = (30 * I)/ (166 * (Y))

        Pc_190 = (6 * I)/ (166 * Y)
        Pt_190 = (30 * I)/ (166 * ((self.n_flanges_bot * self.h_flange + self.h_web + self.n_flanges_top * self.h_flange)-Y))

        Pc_f = min(Pc_166, Pc_190)
        Pt_f = min(Pt_166, Pt_190)
        
        return (Pt_f, Pc_f)

    def calc_P_Shear(self):   
        # Check for Shear both at Glue and Section
        I = self.calc_I()
        Q = self.calc_Q()
        Ps = 6 * I * (self.b_web * self.n_webs)/ Q[1]
        Pg = 2 * I * (self.b_web * self.n_webs)/ Q[0]

        
        return (Pg, Ps)
    
    def calc_P_Plate(self):
        I = self.calc_I()
        Y = self.calc_Centroid() 
        Q = self.calc_Q()

        t = (self.n_flanges_top * self.h_flange)
        h = self.h_web
        a = 530

        # Stress eqn 1
        S_t = (4 * (math.pi**2)* 4000 * (t)**2 ) / ( 12 * 0.96 * (self.b_flange**2)) 

        Pc_166_st = (6 * I)/ (166 * ((self.n_flanges_bot * self.h_flange + self.h_web + self.n_flanges_top * self.h_flange)-Y))
        Pc_190_st = (6 * I)/ (166 * Y)

        Pt_f = min(Pc_166, Pc_190)

        # Stress eqn 3
        S_w = (0.425 * (math.pi**2)* 4000 * (t)**2 ) / ( 12 * 0.96 * (self.b_flange**2))

        Pc_166_sw = (6 * I)/ (166 * ((self.n_flanges_bot * self.h_flange + self.h_web + self.n_flanges_top * self.h_flange)-Y))
        Pc_190_sw = (6 * I)/ (166 * Y)

        Pw_f = min(Pc_166, Pc_190)

        # Shear eqn 4
        Tau =  (  (5 * (math.pi**2)* 4000 * (self.h_flange)**2 )  *  ((t/h)**2 + (t/a)**2)  ) / ( 12 * 0.96)

        PTau_f = Tau * I * (self.b_web * self.n_webs)/ Q[1]
        return(Pt_f, Pw_f, PTau_f)



if __name__ == "__main__":
    n_webs = 2
    n_flanges_top = 1
    n_flanges_bot = 1
    b_web = 1.27
    b_flange = 100 
    h_web = 150
    h_flange = 1.27

    Volume = 1275 * ( h_flange * b_flange * (n_flanges_top + n_flanges_bot) + (n_webs * b_web * h_web) )

    if Volume <= 825000:
        B = Beam(n_webs, n_flanges_top , n_flanges_bot, b_web, b_flange, h_web, h_flange)
        print(B.calc_P_Flex())
        print(B.calc_P_Shear())
        print(B.calc_P_Plate())
        
    else:
        print("Try Again")
