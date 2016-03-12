// ********************************************************************
// Q4 element with bubbles, conforming, 2D
// ********************************************************************

// base function values
static void C_Q_UL4_2D_Funct(double xi, double eta, double *values)
{     
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68;
  double t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81;
  double t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93, t94, t95;
  double t96, t97, t98, t99, t100, t101, t103, t104, t106, t107, t108;
  double t109, t110, t111, t112, t113, t114, t115, t116, t117, t118, t119;
  double t120, t121, t122, t123, t124, t125, t126, t128, t129, t130, t131;
  double t132, t133, t134, t135, t136, t137, t138, t139, t140, t141, t142;
  double t143, t144, t146, t147, t149, t150, t152, t153, t154, t156, t157;
  double t159, t160, t162, t163, t164, t184, t187, t188, t190;

  t1 = 49.0/160.0*xi;
  t2 = xi*xi;
  t3 = t2*t2;
  t4 = t3*xi;
  t5 = 147.0/320.0*t4;
  t6 = 73.0/256.0*t3;
  t7 = t2*xi;
  t8 = 49.0/64.0*t7;
  t9 = 49.0/160.0*eta;
  t10 = eta*eta;
  t11 = t10*t10;
  t12 = 73.0/256.0*t11;
  t13 = t11*eta;
  t14 = 147.0/320.0*t13;
  t15 = t10*eta;
  t16 = 49.0/64.0*t15;
  t17 = t7*eta;
  t18 = 11.0/48.0*t17;
  t19 = 41.0/128.0*t2;
  t20 = 41.0/128.0*t10;
  t21 = xi*eta;
  t22 = 5.0/16.0*t21;
  t23 = t2*t10;
  t24 = 185.0/64.0*t23;
  t25 = -9.0/256.0-t1-t5-t6+t8-t9-t12-t14+t16+t18+t19+t20-t22-t24;
  t26 = t2*eta;
  t27 = 373.0/320.0*t26;
  t28 = t7*t10;
  t29 = 37.0/24.0*t28;
  t30 = t2*t15;
  t31 = 37.0/24.0*t30;
  t32 = t3*t11;
  t33 = 1435.0/768.0*t32;
  t34 = t2*t13;
  t35 = 147.0/320.0*t34;
  t36 = t3*t15;
  t37 = 85.0/192.0*t36;
  t38 = t7*t11;
  t39 = 85.0/192.0*t38;
  t40 = t4*t10;
  t41 = 147.0/320.0*t40;
  t42 = t2*t11;
  t43 = 955.0/384.0*t42;
  t44 = t7*t15;
  t45 = 5.0/48.0*t44;
  t46 = t3*t10;
  t47 = 955.0/384.0*t46;
  t48 = xi*t11;
  t49 = 149.0/192.0*t48;
  t50 = t3*eta;
  t51 = 149.0/192.0*t50;
  t52 = xi*t15;
  t53 = 11.0/48.0*t52;
  t54 = xi*t10;
  t55 = 373.0/320.0*t54;
  t56 = t27-t29-t31-t33+t35+t37+t39+t41+t43+t45+t47-t49-t51+t53+t55;
  t58 = xi/4.0;
  t59 = 9.0/10.0*eta;
  t60 = 21.0/10.0*t34;
  t61 = 35.0/12.0*t48;
  t62 = 5.0*t23;
  t63 = 35.0/6.0*t32;
  t64 = 10.0/3.0*t36;
  t65 = 2.0*t50;
  t66 = 35.0/12.0*t38;
  t67 = 21.0/10.0*t13;
  t68 = 3.0*t15;
  t69 = t7/4.0;
  t70 = -t58-t59+t60-t61-t62-t63+t64-t65+t66-t67+t68+t69;
  t71 = 5.0/2.0*t28;
  t72 = 5.0/2.0*t54;
  t73 = 29.0/10.0*t26;
  t74 = 19.0/3.0*t30;
  t75 = t3/2.0;
  t76 = t2/2.0;
  t77 = 35.0/6.0*t42;
  t78 = 5.0/3.0*t52;
  t79 = 5.0/3.0*t44;
  t80 = 5.0*t46;
  t81 = -t71+t72+t73-t74-t21-t75+t76+t77+t78+t17-t79+t80;
  t83 = 33.0/80.0*eta;
  t84 = 63.0/80.0*t34;
  t85 = 75.0/8.0*t23;
  t86 = 35.0/4.0*t32;
  t87 = 5.0*t36;
  t88 = 3.0*t50;
  t89 = 63.0/80.0*t13;
  t90 = t15/8.0;
  t91 = 273.0/80.0*t26;
  t92 = 41.0/8.0*t30;
  t93 = 3.0/4.0*t3;
  t94 = 15.0/16.0*t2;
  t95 = 15.0/8.0*t10;
  t96 = 175.0/16.0*t42;
  t97 = 35.0/16.0*t11;
  t98 = 15.0/2.0*t46;
  t99 = t83+t84+t85+t86-t87+t88-t89+3.0/16.0-t90-t91+t92+t93-t94-t95-t96+t97-t98;
  t100 = t58-t59+t60+t61-t62-t63+t64-t65-t66-t67+t68-t69;
  t101 = t71-t72+t73-t74+t21-t75+t76+t77-t78-t17+t79+t80;
  t103 = -9.0/256.0+t1+t5-t6-t8-t9-t12-t14+t16-t18+t19+t20+t22-t24;
  t104 = t27+t29-t31-t33+t35+t37-t39-t41+t43-t45+t47+t49-t51-t53-t55;
  t106 = 9.0/10.0*xi;
  t107 = eta/4.0;
  t108 = 2.0*t48;
  t109 = 35.0/12.0*t36;
  t110 = 35.0/12.0*t50;
  t111 = 10.0/3.0*t38;
  t112 = 21.0/10.0*t40;
  t113 = 21.0/10.0*t4;
  t114 = t15/4.0;
  t115 = 3.0*t7;
  t116 = t106-t107+t108-t62-t63+t109-t110-t111-t112+t113+t114-t115;
  t117 = 19.0/3.0*t28;
  t118 = 29.0/10.0*t54;
  t119 = 5.0/2.0*t26;
  t120 = 5.0/2.0*t30;
  t121 = t10/2.0;
  t122 = 5.0*t42;
  t123 = t11/2.0;
  t124 = 5.0/3.0*t17;
  t125 = 35.0/6.0*t46;
  t126 = t117-t118+t119-t120+t21+t121+t122-t123-t52-t124+t79+t125;
  t128 = 33.0/80.0*xi;
  t129 = 3.0*t48;
  t130 = 5.0*t38;
  t131 = 63.0/80.0*t40;
  t132 = 63.0/80.0*t4;
  t133 = t7/8.0;
  t134 = 41.0/8.0*t28;
  t135 = 273.0/80.0*t54;
  t136 = 35.0/16.0*t3;
  t137 = 15.0/8.0*t2;
  t138 = 15.0/16.0*t10;
  t139 = 15.0/2.0*t42;
  t140 = 3.0/4.0*t11;
  t141 = 175.0/16.0*t46;
  t142 = -t128-t129+t85+t86+t130-t131+t132+3.0/16.0+t133-t134+t135+t136-t137-t138-t139+t140-t141;
  t143 = t106+t107+t108-t62-t63-t109+t110-t111-t112+t113-t114-t115;
  t144 = t117-t118-t119+t120-t21+t121+t122-t123+t52+t124-t79+t125;
  t146 = -9.0/256.0+t1+t5-t6-t8+t9-t12+t14-t16+t18+t19+t20-t22-t24;
  t147 = -t27+t29+t31-t33-t35-t37-t39-t41+t43+t45+t47+t49+t51+t53-t55;
  t149 = t58+t59-t60+t61-t62-t63-t64+t65-t66+t67-t68-t69;
  t150 = t71-t72-t73+t74-t21-t75+t76+t77+t78+t17-t79+t80;
  t152 = -t83-t84+t85+t86+t87-t88+t89+3.0/16.0+t90+t91-t92+t93-t94-t95-t96+t97-t98;
  t153 = -t58+t59-t60-t61-t62-t63-t64+t65+t66+t67-t68+t69;
  t154 = -t71+t72-t73+t74+t21-t75+t76+t77-t78-t17+t79+t80;
  t156 = -9.0/256.0-t1-t5-t6+t8+t9-t12+t14-t16-t18+t19+t20+t22-t24;
  t157 = -t27-t29+t31-t33-t35-t37+t39+t41+t43-t45+t47-t49+t51-t53+t55;
  t159 = -t106+t107-t108-t62-t63-t109+t110+t111+t112-t113-t114+t115;
  t160 = -t117+t118-t119+t120+t21+t121+t122-t123-t52-t124+t79+t125;
  t162 = t128+t129+t85+t86-t130+t131-t132+3.0/16.0-t133+t134-t135+t136-t137-t138-t139+t140-t141;
  t163 = -t106-t107-t108-t62-t63+t109-t110+t111+t112-t113+t114+t115;
  t164 = -t117+t118+t119-t120-t21+t121+t122-t123+t52+t124-t79+t125;
  t184 = 525.0/128.0*t10;
  t187 = 525.0/128.0*t2;
  t188 = 1575.0/64.0*t23;
  t190 = 6125.0/256.0*t32;

  values[0] = t25+t56;
  values[1] = t70+t81;
  values[2] = t99;
  values[3] = t100+t101;
  values[4] = t103+t104;
  values[5] = t116+t126;
  values[6] = t142;
  values[7] = t143+t144;
  values[8] = t146+t147;
  values[9] = t149+t150;
  values[10] = t152;
  values[11] = t153+t154;
  values[12] = t156+t157;
  values[13] = t159+t160;
  values[14] = t162;
  values[15] = t163+t164;
  values[16] = -175.0/256.0*t3+75.0/128.0*t2+75.0/128.0*t10+25.0/256.0-175.0/256.0*t11+1225.0/256.0*t32-525.0/128.0*t46-525.0/128.0*t42+225.0/64.0*t23;
  values[17] = -21.0/8.0*eta+735.0/64.0*t15-567.0/64.0*t13-315.0/16.0*t30+693.0/64.0*t26+567.0/64.0*t34-525.0/64.0*t50+525.0/64.0*t36;
  values[18] = -175.0/256.0+t184+1225.0/256.0*t3-3675.0/128.0*t46-t187+t188-875.0/256.0*t11+t190-2625.0/128.0*t42;
  values[19] = -21.0/8.0*xi+693.0/64.0*t54-525.0/64.0*t48+735.0/64.0*t7-315.0/16.0*t28-567.0/64.0*t4+567.0/64.0*t40+525.0/64.0*t38;
  values[20] = 225.0/16.0*t21-225.0/16.0*t52-225.0/16.0*t17+225.0/16.0*t44;
  values[21] = -525.0/64.0*xi+1575.0/32.0*t54-2625.0/64.0*t48+525.0/64.0*t7+2625.0/64.0*t38-1575.0/32.0*t28;
  values[22] = -175.0/256.0+1225.0/256.0*t11-t184+t187-3675.0/128.0*t42+t188-875.0/256.0*t3+t190-2625.0/128.0*t46;
  values[23] = -525.0/64.0*eta-2625.0/64.0*t50+1575.0/32.0*t26+525.0/64.0*t15+2625.0/64.0*t36-1575.0/32.0*t30;
  values[24] = 11025.0/64.0*t23+30625.0/256.0*t32+1225.0/256.0+6125.0/256.0*t3-3675.0/128.0*t2-3675.0/128.0*t10-18375.0/128.0*t42+6125.0/256.0*t11-18375.0/128.0*t46;
  values[25] = -567.0/64.0*eta+945.0/32.0*t15-1323.0/64.0*t13-945.0/32.0*t30+567.0/64.0*t26+1323.0/64.0*t34;
  values[26] = -567.0/64.0*xi+567.0/64.0*t54+945.0/32.0*t7-945.0/32.0*t28-1323.0/64.0*t4+1323.0/64.0*t40;
}


// values of the derivatives in xi direction
static void C_Q_UL4_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
  double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42;
  double t43, t44, t45, t47, t48, t49, t50, t51, t52, t53, t54, t55, t57;
  double t58, t59, t60, t61, t62, t63, t64, t65, t66, t68, t69, t70, t71;
  double t72, t73, t74, t75, t76, t77, t78, t79, t81, t83, t84, t86, t87;
  double t88, t89, t90, t91, t92, t93, t94, t95, t96, t97, t98, t99, t100;
  double t101, t102, t103, t104, t105, t106, t107, t108, t109, t110, t111;
  double t112, t113, t114, t115, t118, t120, t122, t124, t125, t127, t128;
  double t129, t145, t146, t147;

  t1 = eta*eta;
  t2 = t1*t1;
  t4 = xi*t2*eta;
  t5 = 147.0/160.0*t4;
  t6 = 149.0/192.0*t2;
  t7 = xi*t1;
  t8 = 185.0/32.0*t7;
  t9 = xi*xi;
  t10 = t9*xi;
  t11 = t2*t10;
  t12 = 1435.0/192.0*t11;
  t13 = t1*eta;
  t14 = t10*t13;
  t15 = 85.0/48.0*t14;
  t16 = t10*eta;
  t17 = 149.0/48.0*t16;
  t18 = t9*t2;
  t19 = 85.0/64.0*t18;
  t20 = t9*t9;
  t21 = t20*t1;
  t22 = 147.0/64.0*t21;
  t23 = 147.0/64.0*t20;
  t24 = 147.0/64.0*t9;
  t25 = -49.0/160.0+t5-t6-t8-t12+t15-t17+t19+t22-t23+t24;
  t26 = t9*t1;
  t27 = 37.0/8.0*t26;
  t28 = 373.0/320.0*t1;
  t29 = xi*eta;
  t30 = 373.0/160.0*t29;
  t31 = xi*t13;
  t32 = 37.0/12.0*t31;
  t33 = 5.0/16.0*eta;
  t34 = 73.0/64.0*t10;
  t35 = 41.0/64.0*xi;
  t36 = t2*xi;
  t37 = 955.0/192.0*t36;
  t38 = 11.0/48.0*t13;
  t39 = t9*eta;
  t40 = 11.0/16.0*t39;
  t41 = t9*t13;
  t42 = 5.0/16.0*t41;
  t43 = t10*t1;
  t44 = 955.0/96.0*t43;
  t45 = -t27+t28+t30-t32-t33-t34+t35+t37+t38+t40+t42+t44;
  t47 = 21.0/5.0*t4;
  t48 = 35.0/12.0*t2;
  t49 = 10.0*t7;
  t50 = 70.0/3.0*t11;
  t51 = 40.0/3.0*t14;
  t52 = 8.0*t16;
  t53 = 35.0/4.0*t18;
  t54 = 3.0/4.0*t9;
  t55 = 15.0/2.0*t26;
  t57 = 5.0/2.0*t1;
  t58 = 29.0/5.0*t29;
  t59 = 38.0/3.0*t31;
  t60 = 2.0*t10;
  t61 = 35.0/3.0*t36;
  t62 = 5.0/3.0*t13;
  t63 = 3.0*t39;
  t64 = 5.0*t41;
  t65 = 20.0*t43;
  t66 = t57+t58-t59-eta-t60+xi+t61+t62+t63-t64+t65;
  t68 = 63.0/40.0*t4;
  t69 = 75.0/4.0*t7;
  t70 = 35.0*t11;
  t71 = 20.0*t14;
  t72 = 12.0*t16;
  t73 = 273.0/40.0*t29;
  t74 = 41.0/4.0*t31;
  t75 = 3.0*t10;
  t76 = 15.0/8.0*xi;
  t77 = 175.0/8.0*t36;
  t78 = 30.0*t43;
  t79 = t68+t69+t70-t71+t72-t73+t74+t75-t76-t77-t78;
  t81 = -t57+t58-t59+eta-t60+xi+t61-t62-t63+t64+t65;
  t83 = 49.0/160.0+t5+t6-t8-t12+t15-t17-t19-t22+t23-t24;
  t84 = t27-t28+t30-t32+t33-t34+t35+t37-t38-t40-t42+t44;
  t86 = 2.0*t2;
  t87 = 35.0/3.0*t14;
  t88 = 35.0/3.0*t16;
  t89 = 10.0*t18;
  t90 = 21.0/2.0*t21;
  t91 = 21.0/2.0*t20;
  t92 = 9.0*t9;
  t93 = 19.0*t26;
  t94 = 29.0/10.0*t1;
  t95 = 5.0*t29;
  t96 = 5.0*t31;
  t97 = 10.0*t36;
  t98 = 5.0*t39;
  t99 = 70.0/3.0*t43;
  t100 = 9.0/10.0+t86-t49-t50+t87-t88-t89-t90+t91-t92+t93-t94+t95-t96+eta+t97-t13-t98+t64+t99;
  t101 = 3.0*t2;
  t102 = 15.0*t18;
  t103 = 63.0/16.0*t21;
  t104 = 63.0/16.0*t20;
  t105 = 3.0/8.0*t9;
  t106 = 123.0/8.0*t26;
  t107 = 273.0/80.0*t1;
  t108 = 35.0/4.0*t10;
  t109 = 15.0/4.0*xi;
  t110 = 15.0*t36;
  t111 = 175.0/4.0*t43;
  t112 = -33.0/80.0-t101+t69+t70+t102-t103+t104+t105-t106+t107+t108-t109-t110-t111;
  t113 = 9.0/10.0+t86-t49-t50-t87+t88-t89-t90+t91-t92+t93-t94-t95+t96-eta+t97+t13+t98-t64+t99;
  t114 = 49.0/160.0-t5+t6-t8-t12-t15+t17-t19-t22+t23-t24;
  t115 = t27-t28-t30+t32-t33-t34+t35+t37+t38+t40+t42+t44;
  t118 = -t57-t58+t59-eta-t60+xi+t61+t62+t63-t64+t65;
  t120 = -t68+t69+t70+t71-t72+t73-t74+t75-t76-t77-t78;
  t122 = t57-t58+t59+eta-t60+xi+t61-t62-t63+t64+t65;
  t124 = -49.0/160.0-t5-t6-t8-t12-t15+t17+t19+t22-t23+t24;
  t125 = -t27+t28-t30+t32+t33-t34+t35+t37-t38-t40-t42+t44;
  t127 = -9.0/10.0-t86-t49-t50-t87+t88+t89+t90-t91+t92-t93+t94-t95+t96+eta+t97-t13-t98+t64+t99;
  t128 = 33.0/80.0+t101+t69+t70-t102+t103-t104-t105+t106-t107+t108-t109-t110-t111;
  t129 = -9.0/10.0-t86-t49-t50+t87-t88+t89+t90-t91+t92-t93+t94+t95-t96-eta+t97+t13+t98-t64+t99;
  t145 = 525.0/64.0*xi;
  t146 = 1575.0/32.0*t7;
  t147 = 6125.0/64.0*t11;

  values[0] = t25+t45;
  values[1] = -1.0/4.0+t47-t48-t49-t50+t51-t52+t53+t54-t55+t66;
  values[2] = t79;
  values[3] = 1.0/4.0+t47+t48-t49-t50+t51-t52-t53-t54+t55+t81;
  values[4] = t83+t84;
  values[5] = t100;
  values[6] = t112;
  values[7] = t113;
  values[8] = t114+t115;
  values[9] = 1.0/4.0-t47+t48-t49-t50-t51+t52-t53-t54+t55+t118;
  values[10] = t120;
  values[11] = -1.0/4.0-t47-t48-t49-t50-t51+t52+t53+t54-t55+t122;
  values[12] = t124+t125;
  values[13] = t127;
  values[14] = t128;
  values[15] = t129;
  values[16] = -175.0/64.0*t10+75.0/64.0*xi+1225.0/64.0*t11-525.0/32.0*t43-525.0/64.0*t36+225.0/32.0*t7;
  values[17] = -315.0/8.0*t31+693.0/32.0*t29+567.0/32.0*t4-525.0/16.0*t16+525.0/16.0*t14;
  values[18] = 1225.0/64.0*t10-3675.0/32.0*t43-t145+t146+t147-2625.0/64.0*t36;
  values[19] = -21.0/8.0+693.0/64.0*t1-525.0/64.0*t2+2205.0/64.0*t9-945.0/16.0*t26-2835.0/64.0*t20+2835.0/64.0*t21+1575.0/64.0*t18;
  values[20] = 225.0/16.0*eta-225.0/16.0*t13-675.0/16.0*t39+675.0/16.0*t41;
  values[21] = -525.0/64.0+1575.0/32.0*t1-2625.0/64.0*t2+1575.0/64.0*t9+7875.0/64.0*t18-4725.0/32.0*t26;
  values[22] = t145-3675.0/64.0*t36+t146-875.0/64.0*t10+t147-2625.0/32.0*t43;
  values[23] = -2625.0/16.0*t16+1575.0/16.0*t29+2625.0/16.0*t14-1575.0/16.0*t31;
  values[24] = 11025.0/32.0*t7+30625.0/64.0*t11+6125.0/64.0*t10-3675.0/64.0*xi-18375.0/64.0*t36-18375.0/32.0*t43;
  values[25] = -945.0/16.0*t31+567.0/32.0*t29+1323.0/32.0*t4;
  values[26] = -567.0/64.0+567.0/64.0*t1+2835.0/32.0*t9-2835.0/32.0*t26-6615.0/64.0*t20+6615.0/64.0*t21;
}

// values of the derivatives in eta direction
static void C_Q_UL4_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t21, t22, t23, t24, t25, t26, t27, t28, t29;
  double t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42;
  double t43, t44, t45, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56;
  double t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67, t68, t69;
  double t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81, t83;
  double t84, t85, t86, t87, t88, t89, t91, t92, t93, t94, t95, t96, t97;
  double t98, t99, t101, t102, t103, t104, t105, t106, t107, t108, t109;
  double t110, t112, t114, t115, t117, t118, t119, t120, t121, t124, t126;
  double t128, t145, t147, t149;

  t1 = xi*xi;
  t2 = eta*eta;
  t3 = t2*t2;
  t4 = t1*t3;
  t5 = 147.0/64.0*t4;
  t6 = t2*eta;
  t7 = xi*t6;
  t8 = 149.0/48.0*t7;
  t9 = t1*eta;
  t10 = 185.0/32.0*t9;
  t11 = t1*t1;
  t12 = t11*t6;
  t13 = 1435.0/192.0*t12;
  t14 = t2*t11;
  t15 = 85.0/64.0*t14;
  t16 = 149.0/192.0*t11;
  t17 = t1*xi;
  t18 = t17*t6;
  t19 = 85.0/48.0*t18;
  t21 = t11*xi*eta;
  t22 = 147.0/160.0*t21;
  t23 = 147.0/64.0*t3;
  t24 = 147.0/64.0*t2;
  t25 = -49.0/160.0+t5-t8-t10-t13+t15-t16+t19+t22-t23+t24;
  t26 = t17*eta;
  t27 = 37.0/12.0*t26;
  t28 = xi*eta;
  t29 = 373.0/160.0*t28;
  t30 = 373.0/320.0*t1;
  t31 = t1*t2;
  t32 = 37.0/8.0*t31;
  t33 = 5.0/16.0*xi;
  t34 = 41.0/64.0*eta;
  t35 = t1*t6;
  t36 = 955.0/96.0*t35;
  t37 = 73.0/64.0*t6;
  t38 = xi*t2;
  t39 = 11.0/16.0*t38;
  t40 = 11.0/48.0*t17;
  t41 = t17*t2;
  t42 = 5.0/16.0*t41;
  t43 = t11*eta;
  t44 = 955.0/192.0*t43;
  t45 = -t27+t29+t30-t32-t33+t34+t36-t37+t39+t40+t42+t44;
  t47 = 21.0/2.0*t4;
  t48 = 35.0/3.0*t7;
  t49 = 10.0*t9;
  t50 = 70.0/3.0*t12;
  t51 = 10.0*t14;
  t52 = 2.0*t11;
  t53 = 35.0/3.0*t18;
  t54 = 21.0/2.0*t3;
  t55 = 9.0*t2;
  t56 = 5.0*t26;
  t57 = 5.0*t28;
  t58 = 29.0/10.0*t1;
  t59 = 19.0*t31;
  t60 = 70.0/3.0*t35;
  t61 = 5.0*t38;
  t62 = 5.0*t41;
  t63 = 10.0*t43;
  t64 = -9.0/10.0+t47-t48-t49-t50+t51-t52+t53-t54+t55-t56+t57+t58-t59-xi+t60+t61+t17-t62+t63;
  t65 = 63.0/16.0*t4;
  t66 = 75.0/4.0*t9;
  t67 = 35.0*t12;
  t68 = 15.0*t14;
  t69 = 3.0*t11;
  t70 = 63.0/16.0*t3;
  t71 = 3.0/8.0*t2;
  t72 = 273.0/80.0*t1;
  t73 = 123.0/8.0*t31;
  t74 = 15.0/4.0*eta;
  t75 = 175.0/4.0*t35;
  t76 = 35.0/4.0*t6;
  t77 = 15.0*t43;
  t78 = 33.0/80.0+t65+t66+t67-t68+t69-t70-t71-t72+t73-t74-t75+t76-t77;
  t79 = -9.0/10.0+t47+t48-t49-t50+t51-t52-t53-t54+t55+t56-t57+t58-t59+xi+t60-t61-t17+t62+t63;
  t80 = -49.0/160.0+t5+t8-t10-t13+t15-t16-t19-t22-t23+t24;
  t81 = t27-t29+t30-t32+t33+t34+t36-t37-t39-t40-t42+t44;
  t83 = 8.0*t7;
  t84 = 35.0/4.0*t14;
  t85 = 35.0/12.0*t11;
  t86 = 40.0/3.0*t18;
  t87 = 21.0/5.0*t21;
  t88 = 3.0/4.0*t2;
  t89 = 38.0/3.0*t26;
  t91 = 29.0/5.0*t28;
  t92 = 5.0/2.0*t1;
  t93 = 15.0/2.0*t31;
  t94 = 20.0*t35;
  t95 = 2.0*t6;
  t96 = 3.0*t38;
  t97 = 5.0/3.0*t17;
  t98 = 35.0/3.0*t43;
  t99 = -t91+t92-t93+xi+eta+t94-t95-t96-t97+t62+t98;
  t101 = 12.0*t7;
  t102 = 20.0*t18;
  t103 = 63.0/40.0*t21;
  t104 = 41.0/4.0*t26;
  t105 = 273.0/40.0*t28;
  t106 = 15.0/8.0*eta;
  t107 = 30.0*t35;
  t108 = 3.0*t6;
  t109 = 175.0/8.0*t43;
  t110 = -t101+t66+t67+t102-t103-t104+t105-t106-t107+t108-t109;
  t112 = -t91-t92+t93-xi+eta+t94-t95+t96+t97-t62+t98;
  t114 = 49.0/160.0-t5+t8-t10-t13-t15+t16-t19-t22+t23-t24;
  t115 = t27-t29-t30+t32-t33+t34+t36-t37+t39+t40+t42+t44;
  t117 = 9.0/10.0-t47+t48-t49-t50-t51+t52-t53+t54-t55+t56-t57-t58+t59-xi+t60+t61+t17-t62+t63;
  t118 = -33.0/80.0-t65+t66+t67+t68-t69+t70+t71+t72-t73-t74-t75+t76-t77;
  t119 = 9.0/10.0-t47-t48-t49-t50-t51+t52+t53+t54-t55-t56+t57-t58+t59+xi+t60-t61-t17+t62+t63;
  t120 = 49.0/160.0-t5-t8-t10-t13-t15+t16+t19+t22+t23-t24;
  t121 = -t27+t29-t30+t32+t33+t34+t36-t37-t39-t40-t42+t44;
  t124 = t91-t92+t93+xi+eta+t94-t95-t96-t97+t62+t98;
  t126 = t101+t66+t67-t102+t103+t104-t105-t106-t107+t108-t109;
  t128 = t91+t92-t93-xi+eta+t94-t95+t96+t97-t62+t98;
  t145 = 525.0/64.0*eta;
  t147 = 1575.0/32.0*t9;
  t149 = 6125.0/64.0*t12;

  values[0] = t25+t45;
  values[1] = t64;
  values[2] = t78;
  values[3] = t79;
  values[4] = t80+t81;
  values[5] = -1.0/4.0+t83-t49-t50+t84-t85-t86-t87+t88+t89+t99;
  values[6] = t110;
  values[7] = 1.0/4.0+t83-t49-t50-t84+t85-t86-t87-t88+t89+t112;
  values[8] = t114+t115;
  values[9] = t117;
  values[10] = t118;
  values[11] = t119;
  values[12] = t120+t121;
  values[13] = 1.0/4.0-t83-t49-t50-t84+t85+t86+t87-t88-t89+t124;
  values[14] = t126;
  values[15] = -1.0/4.0-t83-t49-t50+t84-t85+t86+t87+t88-t89+t128;
  values[16] = 75.0/64.0*eta-175.0/64.0*t6+1225.0/64.0*t12-525.0/64.0*t43-525.0/32.0*t35+225.0/32.0*t9;
  values[17] = -21.0/8.0+2205.0/64.0*t2-2835.0/64.0*t3-945.0/16.0*t31+693.0/64.0*t1+2835.0/64.0*t4-525.0/64.0*t11+1575.0/64.0*t14;
  values[18] = t145-3675.0/64.0*t43+t147-875.0/64.0*t6+t149-2625.0/32.0*t35;
  values[19] = 693.0/32.0*t28-525.0/16.0*t7-315.0/8.0*t26+567.0/32.0*t21+525.0/16.0*t18;
  values[20] = 225.0/16.0*xi-675.0/16.0*t38-225.0/16.0*t17+675.0/16.0*t41;
  values[21] = 1575.0/16.0*t28-2625.0/16.0*t7+2625.0/16.0*t18-1575.0/16.0*t26;
  values[22] = 1225.0/64.0*t6-t145-3675.0/32.0*t35+t147+t149-2625.0/64.0*t43;
  values[23] = -525.0/64.0-2625.0/64.0*t11+1575.0/32.0*t1+1575.0/64.0*t2+7875.0/64.0*t14-4725.0/32.0*t31;
  values[24] = 11025.0/32.0*t9+30625.0/64.0*t12-3675.0/64.0*eta-18375.0/32.0*t35+6125.0/64.0*t6-18375.0/64.0*t43;
  values[25] = -567.0/64.0+2835.0/32.0*t2-6615.0/64.0*t3-2835.0/32.0*t31+567.0/64.0*t1+6615.0/64.0*t4;
  values[26] = 567.0/32.0*t28-945.0/16.0*t26+1323.0/32.0*t21;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL4_2D_DeriveXiXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67;
  double t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80;
  double t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93;
  double t108, t109;

  t1 = eta*eta;
  t2 = t1*t1;
  t3 = t2*eta;
  t4 = 147.0/160.0*t3;
  t5 = 185.0/32.0*t1;
  t6 = xi*xi;
  t7 = t6*t2;
  t8 = 1435.0/64.0*t7;
  t9 = t1*eta;
  t10 = t6*t9;
  t11 = 85.0/16.0*t10;
  t12 = t6*eta;
  t13 = 149.0/16.0*t12;
  t14 = xi*t2;
  t15 = 85.0/32.0*t14;
  t16 = xi*t6;
  t17 = t16*t1;
  t18 = 147.0/16.0*t17;
  t19 = 147.0/16.0*t16;
  t20 = 147.0/32.0*xi;
  t21 = t1*xi;
  t22 = 37.0/4.0*t21;
  t23 = 373.0/160.0*eta;
  t24 = 37.0/12.0*t9;
  t25 = 219.0/64.0*t6;
  t26 = 955.0/192.0*t2;
  t27 = xi*eta;
  t28 = 11.0/8.0*t27;
  t29 = xi*t9;
  t30 = 5.0/8.0*t29;
  t31 = t1*t6;
  t32 = 955.0/32.0*t31;
  t33 = t4-t5-t8+t11-t13+t15+t18-t19+t20-t22+t23-t24-t25+41.0/64.0+t26+t28+t30+t32;
  t34 = 21.0/5.0*t3;
  t35 = 10.0*t1;
  t36 = 70.0*t7;
  t37 = 40.0*t10;
  t38 = 24.0*t12;
  t39 = 35.0/2.0*t14;
  t40 = 3.0/2.0*xi;
  t41 = 15.0*t21;
  t42 = 29.0/5.0*eta;
  t43 = 38.0/3.0*t9;
  t44 = 6.0*t6;
  t45 = 35.0/3.0*t2;
  t46 = 6.0*t27;
  t47 = 10.0*t29;
  t48 = 60.0*t31;
  t49 = t34-t35-t36+t37-t38+t39+t40-t41+t42-t43-t44+1.0+t45+t46-t47+t48;
  t50 = 63.0/40.0*t3;
  t51 = 75.0/4.0*t1;
  t52 = 105.0*t7;
  t53 = 60.0*t10;
  t54 = 36.0*t12;
  t55 = 273.0/40.0*eta;
  t56 = 41.0/4.0*t9;
  t57 = 9.0*t6;
  t58 = 175.0/8.0*t2;
  t59 = 90.0*t31;
  t60 = t50+t51+t52-t53+t54-t55+t56+t57-15.0/8.0-t58-t59;
  t61 = t34-t35-t36+t37-t38-t39-t40+t41+t42-t43-t44+1.0+t45-t46+t47+t48;
  t62 = t4-t5-t8+t11-t13-t15-t18+t19-t20+t22+t23-t24-t25+41.0/64.0+t26-t28-t30+t32;
  t63 = 35.0*t10;
  t64 = 35.0*t12;
  t65 = 20.0*t14;
  t66 = 42.0*t17;
  t67 = 42.0*t16;
  t68 = 18.0*xi;
  t69 = 38.0*t21;
  t70 = 5.0*eta;
  t71 = 5.0*t9;
  t72 = 10.0*t2;
  t73 = 10.0*t27;
  t74 = 70.0*t31;
  t75 = -t35-t36+t63-t64-t65-t66+t67-t68+t69+t70-t71+t72-t73+t47+t74;
  t76 = 30.0*t14;
  t77 = 63.0/4.0*t17;
  t78 = 63.0/4.0*t16;
  t79 = 3.0/4.0*xi;
  t80 = 123.0/4.0*t21;
  t81 = 105.0/4.0*t6;
  t82 = 15.0*t2;
  t83 = 525.0/4.0*t31;
  t84 = t51+t52+t76-t77+t78+t79-t80+t81-15.0/4.0-t82-t83;
  t85 = -t35-t36-t63+t64-t65-t66+t67-t68+t69-t70+t71+t72+t73-t47+t74;
  t86 = -t4-t5-t8-t11+t13-t15-t18+t19-t20+t22-t23+t24-t25+41.0/64.0+t26+t28+t30+t32;
  t87 = -t34-t35-t36-t37+t38-t39-t40+t41-t42+t43-t44+1.0+t45+t46-t47+t48;
  t88 = -t50+t51+t52+t53-t54+t55-t56+t57-15.0/8.0-t58-t59;
  t89 = -t34-t35-t36-t37+t38+t39+t40-t41-t42+t43-t44+1.0+t45-t46+t47+t48;
  t90 = -t4-t5-t8-t11+t13+t15+t18-t19+t20-t22-t23+t24-t25+41.0/64.0+t26-t28-t30+t32;
  t91 = -t35-t36-t63+t64+t65+t66-t67+t68-t69-t70+t71+t72-t73+t47+t74;
  t92 = t51+t52-t76+t77-t78-t79+t80+t81-15.0/4.0-t82-t83;
  t93 = -t35-t36+t63-t64+t65+t66-t67+t68-t69+t70-t71+t72+t73-t47+t74;
  t108 = 1575.0/32.0*t1;
  t109 = 18375.0/64.0*t7;

  values[0] = t33;
  values[1] = t49;
  values[2] = t60;
  values[3] = t61;
  values[4] = t62;
  values[5] = t75;
  values[6] = t84;
  values[7] = t85;
  values[8] = t86;
  values[9] = t87;
  values[10] = t88;
  values[11] = t89;
  values[12] = t90;
  values[13] = t91;
  values[14] = t92;
  values[15] = t93;
  values[16] = -525.0/64.0*t6+75.0/64.0+3675.0/64.0*t7-1575.0/32.0*t31-525.0/64.0*t2+225.0/32.0*t1;
  values[17] = -315.0/8.0*t9+693.0/32.0*eta+567.0/32.0*t3-1575.0/16.0*t12+1575.0/16.0*t10;
  values[18] = 3675.0/64.0*t6-11025.0/32.0*t31-525.0/64.0+t108+t109-2625.0/64.0*t2;
  values[19] = 2205.0/32.0*xi-945.0/8.0*t21-2835.0/16.0*t16+2835.0/16.0*t17+1575.0/32.0*t14;
  values[20] = -675.0/8.0*t27+675.0/8.0*t29;
  values[21] = 1575.0/32.0*xi+7875.0/32.0*t14-4725.0/16.0*t21;
  values[22] = 525.0/64.0-3675.0/64.0*t2+t108-2625.0/64.0*t6+t109-7875.0/32.0*t31;
  values[23] = -7875.0/16.0*t12+1575.0/16.0*eta+7875.0/16.0*t10-1575.0/16.0*t9;
  values[24] = 11025.0/32.0*t1+91875.0/64.0*t7+18375.0/64.0*t6-3675.0/64.0-18375.0/64.0*t2-55125.0/32.0*t31;
  values[25] = -945.0/16.0*t9+567.0/32.0*eta+1323.0/32.0*t3;
  values[26] = 2835.0/16.0*xi-2835.0/16.0*t21-6615.0/16.0*t16+6615.0/16.0*t17;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL4_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t63, t64, t65, t66, t67, t68;
  double t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80, t81;
  double t82, t83, t84, t85, t87, t88, t89, t91, t92, t93, t95, t108, t109;

  t1 = eta*eta;
  t2 = t1*t1;
  t3 = xi*t2;
  t4 = 147.0/32.0*t3;
  t5 = t1*eta;
  t6 = 149.0/48.0*t5;
  t7 = xi*eta;
  t8 = 185.0/16.0*t7;
  t9 = xi*xi;
  t10 = xi*t9;
  t11 = t10*t5;
  t12 = 1435.0/48.0*t11;
  t13 = t10*t1;
  t14 = 85.0/16.0*t13;
  t15 = 149.0/48.0*t10;
  t16 = t9*t5;
  t17 = 85.0/16.0*t16;
  t18 = t9*t9;
  t19 = t18*eta;
  t20 = 147.0/32.0*t19;
  t21 = t9*eta;
  t22 = 37.0/4.0*t21;
  t23 = 373.0/160.0*eta;
  t24 = 373.0/160.0*xi;
  t25 = t1*xi;
  t26 = 37.0/4.0*t25;
  t27 = xi*t5;
  t28 = 955.0/48.0*t27;
  t29 = 11.0/16.0*t1;
  t30 = 11.0/16.0*t9;
  t31 = t9*t1;
  t32 = 15.0/16.0*t31;
  t33 = t10*eta;
  t34 = 955.0/48.0*t33;
  t35 = t4-t6-t8-t12+t14-t15+t17+t20-t22+t23+t24-t26-5.0/16.0+t28+t29+t30+t32+t34;
  t36 = 21.0*t3;
  t37 = 35.0/3.0*t5;
  t38 = 20.0*t7;
  t39 = 280.0/3.0*t11;
  t40 = 40.0*t13;
  t41 = 8.0*t10;
  t42 = 35.0*t16;
  t43 = 15.0*t21;
  t44 = 5.0*eta;
  t45 = 29.0/5.0*xi;
  t46 = 38.0*t25;
  t47 = 140.0/3.0*t27;
  t48 = 5.0*t1;
  t49 = 3.0*t9;
  t50 = 15.0*t31;
  t51 = 40.0*t33;
  t52 = t36-t37-t38-t39+t40-t41+t42-t43+t44+t45-t46-1.0+t47+t48+t49-t50+t51;
  t53 = 63.0/8.0*t3;
  t54 = 75.0/2.0*t7;
  t55 = 140.0*t11;
  t56 = 60.0*t13;
  t57 = 12.0*t10;
  t58 = 273.0/40.0*xi;
  t59 = 123.0/4.0*t25;
  t60 = 175.0/2.0*t27;
  t61 = 60.0*t33;
  t63 = t36+t37-t38-t39+t40-t41-t42+t43-t44+t45-t46+1.0+t47-t48-t49+t50+t51;
  t64 = t4+t6-t8-t12+t14-t15-t17-t20+t22-t23+t24-t26+5.0/16.0+t28-t29-t30-t32+t34;
  t65 = 8.0*t5;
  t66 = 35.0*t13;
  t67 = 35.0/3.0*t10;
  t68 = 40.0*t16;
  t69 = 21.0*t19;
  t70 = 38.0*t21;
  t71 = 29.0/5.0*eta;
  t72 = 5.0*xi;
  t73 = 15.0*t25;
  t74 = 40.0*t27;
  t75 = 3.0*t1;
  t76 = 5.0*t9;
  t77 = 140.0/3.0*t33;
  t78 = t65-t38-t39+t66-t67-t68-t69+t70-t71+t72-t73+1.0+t74-t75-t76+t50+t77;
  t79 = 12.0*t5;
  t80 = 60.0*t16;
  t81 = 63.0/8.0*t19;
  t82 = 123.0/4.0*t21;
  t83 = 273.0/40.0*eta;
  t84 = 60.0*t27;
  t85 = 175.0/2.0*t33;
  t87 = t65-t38-t39-t66+t67-t68-t69+t70-t71-t72+t73-1.0+t74+t75+t76-t50+t77;
  t88 = -t4+t6-t8-t12-t14+t15-t17-t20+t22-t23-t24+t26-5.0/16.0+t28+t29+t30+t32+t34;
  t89 = -t36+t37-t38-t39-t40+t41-t42+t43-t44-t45+t46-1.0+t47+t48+t49-t50+t51;
  t91 = -t36-t37-t38-t39-t40+t41+t42-t43+t44-t45+t46+1.0+t47-t48-t49+t50+t51;
  t92 = -t4-t6-t8-t12-t14+t15+t17+t20-t22+t23-t24+t26+5.0/16.0+t28-t29-t30-t32+t34;
  t93 = -t65-t38-t39-t66+t67+t68+t69-t70+t71-t72+t73+1.0+t74-t75-t76+t50+t77;
  t95 = -t65-t38-t39+t66-t67+t68+t69-t70+t71+t72-t73-1.0+t74+t75+t76-t50+t77;
  t108 = 1575.0/16.0*t7;
  t109 = 6125.0/16.0*t11;

  values[0] = t35;
  values[1] = t52;
  values[2] = t53+t54+t55-t56+t57-t58+t59-t60-t61;
  values[3] = t63;
  values[4] = t64;
  values[5] = t78;
  values[6] = -t79+t54+t55+t80-t81-t82+t83-t84-t85;
  values[7] = t87;
  values[8] = t88;
  values[9] = t89;
  values[10] = -t53+t54+t55+t56-t57+t58-t59-t60-t61;
  values[11] = t91;
  values[12] = t92;
  values[13] = t93;
  values[14] = t79+t54+t55-t80+t81+t82-t83-t84-t85;
  values[15] = t95;
  values[16] = 1225.0/16.0*t11-525.0/16.0*t33-525.0/16.0*t27+225.0/16.0*t7;
  values[17] = -945.0/8.0*t25+693.0/32.0*xi+2835.0/32.0*t3-525.0/16.0*t10+1575.0/16.0*t13;
  values[18] = -3675.0/16.0*t33+t108+t109-2625.0/16.0*t27;
  values[19] = 693.0/32.0*eta-525.0/16.0*t5-945.0/8.0*t21+2835.0/32.0*t19+1575.0/16.0*t16;
  values[20] = 225.0/16.0-675.0/16.0*t1-675.0/16.0*t9+2025.0/16.0*t31;
  values[21] = 1575.0/16.0*eta-2625.0/16.0*t5+7875.0/16.0*t16-4725.0/16.0*t21;
  values[22] = -3675.0/16.0*t27+t108+t109-2625.0/16.0*t33;
  values[23] = -2625.0/16.0*t10+1575.0/16.0*xi+7875.0/16.0*t13-4725.0/16.0*t25;
  values[24] = 11025.0/16.0*t7+30625.0/16.0*t11-18375.0/16.0*t27-18375.0/16.0*t33;
  values[25] = -2835.0/16.0*t25+567.0/32.0*xi+6615.0/32.0*t3;
  values[26] = 567.0/32.0*eta-2835.0/16.0*t21+6615.0/32.0*t19;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL4_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
  double t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54;
  double t55, t56, t57, t58, t59, t60, t61, t62, t63, t64, t65, t66, t67;
  double t68, t69, t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80;
  double t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92, t93;
  double t107, t109;

  t1 = xi*xi;
  t2 = eta*eta;
  t3 = t2*eta;
  t4 = t1*t3;
  t5 = 147.0/16.0*t4;
  t6 = xi*t2;
  t7 = 149.0/16.0*t6;
  t8 = 185.0/32.0*t1;
  t9 = t1*t1;
  t10 = t9*t2;
  t11 = 1435.0/64.0*t10;
  t12 = t9*eta;
  t13 = 85.0/32.0*t12;
  t14 = t1*xi;
  t15 = t14*t2;
  t16 = 85.0/16.0*t15;
  t17 = xi*t9;
  t18 = 147.0/160.0*t17;
  t19 = 147.0/16.0*t3;
  t20 = 147.0/32.0*eta;
  t21 = 37.0/12.0*t14;
  t22 = 373.0/160.0*xi;
  t23 = t1*eta;
  t24 = 37.0/4.0*t23;
  t25 = t1*t2;
  t26 = 955.0/32.0*t25;
  t27 = 219.0/64.0*t2;
  t28 = xi*eta;
  t29 = 11.0/8.0*t28;
  t30 = t14*eta;
  t31 = 5.0/8.0*t30;
  t32 = 955.0/192.0*t9;
  t33 = t5-t7-t8-t11+t13+t16+t18-t19+t20-t21+t22-t24+41.0/64.0+t26-t27+t29+t31+t32;
  t34 = 42.0*t4;
  t35 = 35.0*t6;
  t36 = 10.0*t1;
  t37 = 70.0*t10;
  t38 = 20.0*t12;
  t39 = 35.0*t15;
  t40 = 42.0*t3;
  t41 = 18.0*eta;
  t42 = 5.0*t14;
  t43 = 5.0*xi;
  t44 = 38.0*t23;
  t45 = 70.0*t25;
  t46 = 10.0*t28;
  t47 = 10.0*t30;
  t48 = 10.0*t9;
  t49 = t34-t35-t36-t37+t38+t39-t40+t41-t42+t43-t44+t45+t46-t47+t48;
  t50 = 63.0/4.0*t4;
  t51 = 75.0/4.0*t1;
  t52 = 105.0*t10;
  t53 = 30.0*t12;
  t54 = 63.0/4.0*t3;
  t55 = 3.0/4.0*eta;
  t56 = 123.0/4.0*t23;
  t57 = 525.0/4.0*t25;
  t58 = 105.0/4.0*t2;
  t59 = 15.0*t9;
  t60 = t50+t51+t52-t53-t54-t55+t56-15.0/4.0-t57+t58-t59;
  t61 = t34+t35-t36-t37+t38-t39-t40+t41+t42-t43-t44+t45-t46+t47+t48;
  t62 = t5+t7-t8-t11+t13-t16-t18-t19+t20+t21-t22-t24+41.0/64.0+t26-t27-t29-t31+t32;
  t63 = 24.0*t6;
  t64 = 35.0/2.0*t12;
  t65 = 40.0*t15;
  t66 = 21.0/5.0*t17;
  t67 = 3.0/2.0*eta;
  t68 = 38.0/3.0*t14;
  t69 = 29.0/5.0*xi;
  t70 = 15.0*t23;
  t71 = 60.0*t25;
  t72 = 6.0*t2;
  t73 = 6.0*t28;
  t74 = 35.0/3.0*t9;
  t75 = t63-t36-t37+t64-t65-t66+t67+t68-t69-t70+1.0+t71-t72-t73+t47+t74;
  t76 = 36.0*t6;
  t77 = 60.0*t15;
  t78 = 63.0/40.0*t17;
  t79 = 41.0/4.0*t14;
  t80 = 273.0/40.0*xi;
  t81 = 90.0*t25;
  t82 = 9.0*t2;
  t83 = 175.0/8.0*t9;
  t84 = -t76+t51+t52+t77-t78-t79+t80-15.0/8.0-t81+t82-t83;
  t85 = t63-t36-t37-t64-t65-t66-t67+t68-t69+t70+1.0+t71-t72+t73-t47+t74;
  t86 = -t5+t7-t8-t11-t13-t16-t18+t19-t20+t21-t22+t24+41.0/64.0+t26-t27+t29+t31+t32;
  t87 = -t34+t35-t36-t37-t38-t39+t40-t41+t42-t43+t44+t45+t46-t47+t48;
  t88 = -t50+t51+t52+t53+t54+t55-t56-15.0/4.0-t57+t58-t59;
  t89 = -t34-t35-t36-t37-t38+t39+t40-t41-t42+t43+t44+t45-t46+t47+t48;
  t90 = -t5-t7-t8-t11-t13+t16+t18+t19-t20-t21+t22+t24+41.0/64.0+t26-t27-t29-t31+t32;
  t91 = -t63-t36-t37-t64+t65+t66-t67-t68+t69+t70+1.0+t71-t72-t73+t47+t74;
  t92 = t76+t51+t52-t77+t78+t79-t80-15.0/8.0-t81+t82-t83;
  t93 = -t63-t36-t37+t64+t65+t66+t67-t68+t69-t70+1.0+t71-t72+t73-t47+t74;
  t107 = 1575.0/32.0*t1;
  t109 = 18375.0/64.0*t10;

  values[0] = t33;
  values[1] = t49;
  values[2] = t60;
  values[3] = t61;
  values[4] = t62;
  values[5] = t75;
  values[6] = t84;
  values[7] = t85;
  values[8] = t86;
  values[9] = t87;
  values[10] = t88;
  values[11] = t89;
  values[12] = t90;
  values[13] = t91;
  values[14] = t92;
  values[15] = t93;
  values[16] = 75.0/64.0-525.0/64.0*t2+3675.0/64.0*t10-525.0/64.0*t9-1575.0/32.0*t25+225.0/32.0*t1;
  values[17] = 2205.0/32.0*eta-2835.0/16.0*t3-945.0/8.0*t23+2835.0/16.0*t4+1575.0/32.0*t12;
  values[18] = 525.0/64.0-3675.0/64.0*t9+t107-2625.0/64.0*t2+t109-7875.0/32.0*t25;
  values[19] = 693.0/32.0*xi-1575.0/16.0*t6-315.0/8.0*t14+567.0/32.0*t17+1575.0/16.0*t15;
  values[20] = -675.0/8.0*t28+675.0/8.0*t30;
  values[21] = 1575.0/16.0*xi-7875.0/16.0*t6+7875.0/16.0*t15-1575.0/16.0*t14;
  values[22] = 3675.0/64.0*t2-525.0/64.0-11025.0/32.0*t25+t107+t109-2625.0/64.0*t9;
  values[23] = 1575.0/32.0*eta+7875.0/32.0*t12-4725.0/16.0*t23;
  values[24] = 11025.0/32.0*t1+91875.0/64.0*t10-3675.0/64.0-55125.0/32.0*t25+18375.0/64.0*t2-18375.0/64.0*t9;
  values[25] = 2835.0/16.0*eta-6615.0/16.0*t3-2835.0/16.0*t23+6615.0/16.0*t4;
  values[26] = 567.0/32.0*xi-945.0/16.0*t14+1323.0/32.0*t17;
}


// ***********************************************************************

TBaseFunct2D *BF_C_Q_UL4_2D_Obj = new TBaseFunct2D
        (27, BF_C_Q_UL4_2D, BFUnitSquare, 
         C_Q_UL4_2D_Funct, C_Q_UL4_2D_DeriveXi,
         C_Q_UL4_2D_DeriveEta, C_Q_UL4_2D_DeriveXiXi,
         C_Q_UL4_2D_DeriveXiEta, C_Q_UL4_2D_DeriveEtaEta, 5, 4,
         0, NULL);
