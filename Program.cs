using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PullingPoZavodu3
{
    internal class Program
    {
        static void Main(string[] args)
        {
        }
    }

    public class Functions
    {
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, толщЛенты, верхДопТолщ, нижДопТолщ
        //На выход: d0КолпНешт
        public static double d0KolpNesht(double outerDiameter, double outerDiameterTolerance, double innerDiameter, double innerDiameterTolerance, double height, double tapeThickness, double topAdditionalThickness, double bottomAdditionalThickness)
        {
            double dn = outerDiameter + outerDiameterTolerance / 2;
            double dw = innerDiameter + innerDiameterTolerance / 2;
            double h = height;
            double s = tapeThickness + (topAdditionalThickness + bottomAdditionalThickness) / 2;
            double pi = Math.PI;
            double v = pi * (dn * dn - dw * dw) * (h - s) / 4 + pi * dw * dw * s / 4;
            double d0KolpNesht = Math.Sqrt(4 * v / pi / s);
            d0KolpNesht = Math.Round(d0KolpNesht, 2);
            return d0KolpNesht;
        }

        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, толщЛенты, верхДопТолщ, нижДопТолщ, радиус
        //На выход: d0КолпШт
        public static double d0KolpSht(double outerDiameter, double outerDiameterTolerance, double innerDiameter, double innerDiameterTolerance, double height, double tapeThickness, double topAdditionalThickness, double bottomAdditionalThickness, double radius)
        {
            double dn = outerDiameter + outerDiameterTolerance / 2;
            double dw = innerDiameter + innerDiameterTolerance / 2;
            double h = height;
            double s = tapeThickness + (topAdditionalThickness + bottomAdditionalThickness) / 2;
            double r = radius;
            double pi = Math.PI;
            double vcil = pi * (dn * dn - dw * dw) * (h - s - r) / 4;
            double vdno = pi * s * (dw - 2 * r) * (dw - 2 * r) / 4;
            double vtor = pi * pi * (dw / 2 - r) * ((r + s) * (r + s) - r * r) / 2;
            double v = vcil + vdno + vtor;
            double d0KolpSht = Math.Sqrt(4 * v / pi / s);
            d0KolpSht = Math.Round(d0KolpSht, 2);
            return d0KolpSht;
        }

        //На вход: md0, d0, ms0, sl, delsw, delsn, dn, dw, deldn, deldw, s0, m1
        //На выход: md, ms, dsr, s, dm, dk
        public static void KoeffWyt()
        {
            double md0 = float.Parse(Console.ReadLine());
            double d0 = float.Parse(Console.ReadLine());
            double ms0 = float.Parse(Console.ReadLine());
            double sl = float.Parse(Console.ReadLine());
            double delsw = float.Parse(Console.ReadLine());
            double delsn = float.Parse(Console.ReadLine());
            double dn = float.Parse(Console.ReadLine());
            double dw = float.Parse(Console.ReadLine());
            double deldn = float.Parse(Console.ReadLine());
            double deldw = float.Parse(Console.ReadLine());
            double s0 = sl + delsw / 2 + delsn / 2;

            // Calculate the initial values
            Console.WriteLine("Введи  кофф. свертки: ");
            double m1 = float.Parse(Console.ReadLine());
            double mu1 = 1;
            int n = 0;
            if (md0 - m1 >= 0)
            {
                m1 = md0;
                mu1 = ms0;
                n = 1;
            }
            else
            {
                Console.WriteLine("Введи  кофф.утонен. на свертке: ");
                mu1 = float.Parse(Console.ReadLine());

                Console.WriteLine("Введи средний кофф. вытяжки по операциям: ");
                double mdsr = float.Parse(Console.ReadLine());

                double nsr = (Math.Log(md0) - Math.Log(m1)) / Math.Log(mdsr) + 1;
                n = (int)Math.Round(nsr);
                if ((nsr - n) > 0.3)
                {
                    n = n + 1;
                }
            }
            Console.WriteLine("Количество вытяжек: " + (n - 1));
            double[] md = new double[n];
            double[] dsr = new double[n];
            double[] ms = new double[n];
            double[] s = new double[n];
            double[] dm = new double[n];
            double[] dp = new double[n];
            bool isAverageAttenuationCoefficient = false;
            double x, y;
            if (n > 1)
            {
                Console.WriteLine("Коэф.утон.средн.по вытяж?");
                string vbYesNo = Console.ReadLine();
                if (vbYesNo == "Yes")
                {
                    isAverageAttenuationCoefficient = true;
                }
                else
                {
                    isAverageAttenuationCoefficient = false;
                }
                double mssred = Math.Pow(ms0 / mu1, 1 / (n - 1));
            }

            for (int k = 1; k <= n; k++)
            {
                if (k == 1)
                {
                    md[k] = m1;
                    ms[k] = mu1;
                }
                else if (k < n && isAverageAttenuationCoefficient == false)
                {
                    Console.WriteLine("Введи коэфф. вытяжки: " + (k - 1));
                    md[k] = Convert.ToDouble(Console.ReadLine());
                    Console.WriteLine("Введи коэфф. утонен:" + (k - 1));
                    ms[k] = Convert.ToDouble(Console.ReadLine());
                }
                else if (k < n && isAverageAttenuationCoefficient == true)
                {
                    Console.WriteLine("Введи коэфф. вытяжки: " + (k - 1));
                    md[k] = Convert.ToDouble(Console.ReadLine());
                    //видимо так?????
                    double mssred = Math.Pow(ms0 / mu1, 1 / (n - 1));
                    ms[k] = mssred;
                }
                else if (k == n && n > 1)
                {
                    md[k] = md0 / x;
                    ms[k] = ms0 / y;
                }

                if (k == 1)
                {
                    x = md[k];
                    dsr[k] = d0 * md[k];
                    y = ms[k];
                    s[k] = s0 * ms[k];
                }
                else
                {
                    x = x * md[k];
                    dsr[k] = dsr[k - 1] * md[k];
                    y = y * ms[k];
                    s[k] = s[k - 1] * ms[k];
                }

                if ((k == 1) && (ms[k] == 1) && n > 1)
                {
                    dm[k] = dsr[k] + s[k];
                    dp[k] = dsr[k] - s0 / Math.Sqrt(md[k]);
                }
                else if (k == 1 && n == 1)
                {
                    dm[k] = dn + deldn / 2;
                    dp[k] = dw + deldw / 2;
                }
                else
                {
                    dm[k] = dsr[k] + s[k];
                    dp[k] = dsr[k] - s[k];
                }

                if ((k == n) && (ms[k] == 1) && n > 1)
                {
                    dm[k] = dn + deldn / 2;
                    dp[k] = dw + deldw / 2;
                    dsr[k] = (dm[k] + dp[k]) / 2;
                    md[k] = dsr[k] / dsr[k - 1];
                }

                dm[k] = Math.Round(dm[k], 2);
                dp[k] = Math.Round(dp[k], 2);
                s[k] = Math.Round(s[k], 2);

                Console.WriteLine("  md = " + md[k]);
                Console.WriteLine("  dsr = " + dsr[k]);
                Console.WriteLine("  ms = " + ms[k]);
                Console.WriteLine("  s = " + s[k]);
                Console.WriteLine("  dm = " + dm[k]);
                Console.WriteLine("  dp = " + dp[k]);
            }
        }

        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ,  
        //На выход: 
        public static void DiamKru()
        {
            double outerDiameter = double.Parse(Console.ReadLine());
            double outerDiameterTolerance = double.Parse(Console.ReadLine());
            double innerDiameter = double.Parse(Console.ReadLine());
            double innerDiameterTolerance = double.Parse(Console.ReadLine());
            double height = double.Parse(Console.ReadLine());
            double radius = double.Parse(Console.ReadLine());
            double tapeThickness = double.Parse(Console.ReadLine());
            double topAdditionalThickness = double.Parse(Console.ReadLine());
            double bottomAdditionalThickness = double.Parse(Console.ReadLine());

            Console.Write("Дно штампуется? - ");
            string type = Console.ReadLine();
            double dz;
            if (type == "Yes")
            {
                dz = d0KolpSht(outerDiameter, outerDiameterTolerance, innerDiameter, innerDiameterTolerance, height, tapeThickness, topAdditionalThickness, bottomAdditionalThickness, radius);
            }
            else
            {
                dz = d0KolpNesht(outerDiameter, outerDiameterTolerance, innerDiameter, innerDiameterTolerance, height, tapeThickness, topAdditionalThickness, bottomAdditionalThickness);

            }

            Console.WriteLine("Устраивает диаметр кружка dz = " + dz + " ?");
            string answear = Console.ReadLine();
            if (answear == "Yes")
            {
                Console.WriteLine("диаметр кружка dz = " + dz);
            }
            else
            {
                Console.Write("Введи диаметр кружка сам: ");
                dz = Convert.ToDouble(Console.ReadLine());
            }

            double md0 = (outerDiameter + outerDiameterTolerance / 2 + innerDiameter + innerDiameterTolerance / 2) / 2 / dz;
            double sd = tapeThickness / dz;
            double ms = ((outerDiameter + outerDiameterTolerance / 2 - (innerDiameter + innerDiameterTolerance / 2)) / 2) / (tapeThickness + topAdditionalThickness / 2 + bottomAdditionalThickness / 2);
            if (ms > 1)
            {
                ms = 1;
            }
            Console.WriteLine("md0 = " + md0);
            Console.WriteLine("sd = " + sd);
            Console.WriteLine("ms = " + ms);

            //и тут они блять решили их переименовать нахуй!
            double dn = outerDiameter;
            double deln = outerDiameterTolerance;
            double dw = innerDiameter;
            double delw = innerDiameterTolerance;
            double r = radius;
            double t = tapeThickness;
            double tw = topAdditionalThickness;
            double tn = bottomAdditionalThickness;

            double dp = 0;
            if (type == "Yes")
            {
                double dm = dn + deln / 2;
                dp = dw + delw / 2;
                double s = t + tw;
                double hmax = hh(dz, s, dp, dm, r);
                hmax = Math.Round(hmax, 2);
                Console.WriteLine("Высота детали");
                Console.WriteLine("Инструм. неизнош.");
                Console.WriteLine("Hmax=" + hmax);
                s = t + tn;
                double hmin = hh(dz, s, dp, dm, r);
                hmin = Math.Round(hmin, 2);
                Console.WriteLine("Hmin=" + hmin);
                Console.WriteLine("Инструм.изнош.");
                dm = dn;
                dp = dw;
                s = t + tw;
                hmax = hh(dz, s, dp, dm, r);
                hmax = Math.Round(hmax, 2);
                Console.WriteLine("Hmax=" + hmax);
                s = t + tn;
                hmin = hh(dz, s, dp, dm, r);
                hmin = Math.Round(hmin, 2);
                Console.WriteLine("Hmin=" + hmin);
            }
            else
            {
                double dm = dn + deln / 2;
                dp = dw + delw / 2;
                double s = t + tw;
                double hmax = hn(dz, s, dp, dm);
                hmax = Math.Round(hmax, 2);
                s = t + tn;
                double hmin = hn(dz, s, dp, dm);
                hmin = Math.Round(hmin, 2);
                Console.WriteLine("Высота детали");
                Console.WriteLine("Инструм. неизнош.");
                Console.WriteLine("Hmax=" + hmax);
                Console.WriteLine("Hmin=" + hmin);
                Console.WriteLine("Инструм.изнош.");
                dm = dn;
                dp = dw;
                s = t + tw;
                hmax = hn(dz, s, dp, dm);
                hmax = Math.Round(hmax, 2);
                hmin = hmax / s * (t + tn);
                hmin = Math.Round(hmin, 2);
                Console.WriteLine("Hmax=" + hmax);
                Console.WriteLine("Hmin=" + hmin);
            }

            if (answear == "No" && type == "Yes")
            {
                double dm = dn + deln / 2;
                dp = dw + delw / 2;
                double s = t + (tw + tn) / 2;
                double hnew = hh(dz, s, dp, dm, r);
                hnew = Math.Round(hnew, 2);
                Console.WriteLine("Hnew=" + hnew);
            }
            else if (answear == "No" && type == "No")
            {
                double dm = dn + deln / 2;
                dp = dw + delw / 2;
                double s = t + (tw + tn) / 2;
                double hnew = hn(dz, s, dp, dm);
                hnew = Math.Round(hnew, 2);
                Console.WriteLine("Hnew=" + hnew);
            }

            if ((dz - dp) / 2 < t)
            {
                Console.WriteLine("штамп. сдвигом");
                ВырВытяжкаСдвигом();
            }
        }

        //На вход: диаметр матр. предыдущей вытяжки, диаметр пуанс. предыдущ. вытяжки, диам. пуанс. вытяжки для которой расчет, степень деформац. на данной операции 
        //На выход: диам. матриц(кроме 1-й и окончат.) по отчету ЛВМИ 1960
        public static double DmSmirAl(double dmf, double dpf, double dps, double stdef)
        {
            double a = Math.Sqrt(3) / 2;
            double b = 1 - 1 / Math.Sqrt(3);
            double c = a * (stdef - b * Math.Log(dpf / dps));
            double d = Math.Exp(c);
            return Math.Round(Math.Sqrt(dps * dps + (dmf * dmf - dpf * dpf) / d), 2);
        }

        //На вход: диам. кружка, исходн. толщ. ленты, нар. диам. изделия, внутр. диам изделия
        //На выход: суммар. степень деформ по Смирнову-Аляеву
        public static double StepDefSum(double d0, double s0, double dn, double dw)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(d0 / dw) + b * Math.Log(4 * d0 * s0 / (dn * dn - dw * dw));
        }

        //На вход: диам. кружка, исходн. толщ. ленты, диам. матр. на свертке, диам. пуанс. на свертке
        //На выход: степ. деформ. на свертке по Смирнову-Аляеву
        public static double StepDefSwert(double d0, double s0, double dn, double dw)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(d0 / dw) + b * Math.Log(4 * d0 * s0 / (dn * dn - dw * dw));
        }

        //На вход: диам. матр предыдущ. вытяжки, диам. матр. текущей вытяжки, дам. пуанс. предыдущ. вытяж, дам. пуанс. текущей вытяжки
        //На выход: степ. деформ. на вытяжке по Смирнову-Аляеву
        public static double StepDefWyt(double dmf, double dms, double dpf, double dps)
        {
            double a = 1 - 1 / Math.Sqrt(3);
            double b = 2 / Math.Sqrt(3);
            return a * Math.Log(dpf / dps) + b * Math.Log((dmf * dmf - dpf * dpf) / (dms * dms - dps * dps));
        }

        //На вход:
        //На выход:
        public static double hh(double dz, double s, double dp, double dm, double r)
        {
            double Pi = Math.PI;
            return (dz * dz * s - (dp - 2 * r) * s * (Pi * (2 * r + s) + dp - 2 * r)) / (dm * dm - dp * dp) + s + r;
        }

        //На вход:
        //На выход:
        public static double hn(double dz, double s, double dp, double dm)
        {
            return s * ((dz * dz - dp * dp) / (dm * dm - dp * dp) + 1);
        }
    }


    public class SmirAl
    {
        //хз как правильно надо заполнять эти массивы!
        public static void RadPuans()
        {
            double r = double.Parse(Console.ReadLine());
            int n = int.Parse(Console.ReadLine());

            double[] md = new double[n];
            double[] dsr = new double[n];
            double[] rp = new double[n];

            for (int i = 0; i < n; i++)
            {
                dsr[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
            }

            for (int k = 1; k < n; k++)
            {
                if (k == n)
                {
                    rp[k] = r;
                }
                else
                {
                    rp[k] = dsr[k] * (1 - md[k + 1]) / 2;
                    rp[k] = Math.Round(rp[k], 1);
                }
                Console.WriteLine(rp[k]);
            }
        }

        // непонятно с циклами в массивах!
        public static void WysotaPoluf()
        {
            double r = double.Parse(Console.ReadLine());
            int n = int.Parse(Console.ReadLine());
            double[] dsr = new double[n];
            double[] md = new double[n];
            double[] rp = new double[n];
            double[] ms = new double[n];
            double[] h = new double[n];
            double[] x = new double[n];
            double[] y = new double[n];

            for (int i = 0; i < n; i++)
            {
                dsr[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
                ms[i] = double.Parse(Console.ReadLine());
                rp[i] = double.Parse(Console.ReadLine());
            }

            for (int k = 1; k < n; k++)
            {
                x[0] = 1;
                y[0] = 1;
                x[k] = x[k - 1] * ms[k];
                y[k] = y[k - 1] * md[k];
                h[k] = 0.25 * dsr[k] * (1 / y[k] * y[k] - 1 - 2.28 * rp[k] / dsr[k] + 0.56 * (rp[k] / dsr[k]) * (rp[k] / dsr[k])) / x[k] + rp[k];
                h[k] = Math.Round(h[k], 1);
                Console.WriteLine(h[k]);
            }
        }

        public static void ВысотаОткуса()
        {
            double dmp = double.Parse(Console.ReadLine());
            double dpp = double.Parse(Console.ReadLine());
            double dmo = double.Parse(Console.ReadLine());
            double dpo = double.Parse(Console.ReadLine());
            double al = double.Parse(Console.ReadLine());
            double d0 = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double Pi = Math.PI;
            double h = double.Parse(Console.ReadLine());

            double v0 = Pi * d0 * d0 * s / 4;
            double vd = Pi * ((dmo * dmo - dpo * dpo) * h + dpo * dpo * s) / 4;
            double alf = Pi * al / 180;
            double hk = (dmp - dmo) / (2 * Math.Tan(alf));
            hk = Math.Round(hk, 1);
            double lk = (dmp - dmo) / (2 * Math.Sin(alf));
            lk = Math.Round(lk, 1);
            double sp = (dmp - dpp) / 2;
            double vkon = Pi * lk * (dmp - sp) * (dmp - dpp) / 4 + Pi * lk * (dmo - sp) * (dmp - dpp) / 4;
            double votk = v0 - vd;
            if (votk < 0)
            {
                Console.WriteLine("нет металла на откус");
                return;
            }
            if (vkon >= votk)
            {
                double ds = dmo - sp;
                double a = Pi * Math.Sin(alf) * sp;
                double b = Pi * ds * sp;
                double c = -votk;
                if (b * b - 4 * a * c < 0)
                {
                    Console.WriteLine("наверно неправильно введены данные");
                    return;
                }
                double d = Math.Sqrt(b * b - 4 * a * c);
                double l = (d - b) / (2 * a);
                l = Math.Round(l, 1);
                Console.WriteLine("длина откуса=" + l);
            }
            else
            {
                double vc = votk - vkon;
                double hc = 4 * vc / (Pi * (dmp * dmp - dpp * dpp));
                hc = Math.Round(hc, 1);
                Console.WriteLine("высота цил.части откус=" + hc);
                Console.WriteLine("высот. конич.части отк=" + hk);
                double hs = hk + hc;
                Console.WriteLine("суммар высота откуса=" + hs);
            }
        }

        // непонятно с циклами в массивах!
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ
        //На выход: 
        public static void КоэфИзвестнИнструм()
        {
            int n = int.Parse(Console.ReadLine());
            double d0 = double.Parse(Console.ReadLine());
            double dn = double.Parse(Console.ReadLine());
            double deln = double.Parse(Console.ReadLine());
            double dw = double.Parse(Console.ReadLine());
            double delw = double.Parse(Console.ReadLine());
            double r = double.Parse(Console.ReadLine());
            double t = double.Parse(Console.ReadLine());
            double tw = double.Parse(Console.ReadLine());
            double tn = double.Parse(Console.ReadLine());

            double[] dsr = new double[n];
            double[] md = new double[n];
            double[] ms = new double[n];
            double[] s = new double[n];
            double[] dm = new double[n];
            double[] dp = new double[n];
            double[] rp = new double[n];
            double[] h = new double[n];
            double[] x = new double[n];
            double[] y = new double[n];

            for (int i = 1; i < n; i++)
            {
                rp[i] = double.Parse(Console.ReadLine());
                dm[i] = double.Parse(Console.ReadLine());
                dp[i] = double.Parse(Console.ReadLine());
            }

            for (int i = 1; i < n; i++)
            {
                s[i] = (dm[i] - dp[i]) / 2;
                if (i == 1 && (dm[i] - dp[i]) / 2 > t)
                {
                    s[i] = (t + tw / 2 + tn / 2);
                }
                dsr[i] = (dm[i] + dp[i]) / 2;
                if (i == 1)
                {
                    md[i] = dsr[i] / d0;
                    ms[i] = s[i] / (t + tw / 2 + tn / 2);
                }
                else
                {
                    md[i] = dsr[i] / dsr[i - 1];
                    ms[i] = s[i] / s[i - 1];
                }
                if (ms[i] > 1)
                {
                    ms[i] = 1;
                }
                Console.WriteLine(md[i]);
                Console.WriteLine(ms[i]);
                Console.WriteLine(dsr[i]);
                Console.WriteLine(s[i]);
            }
            double md0 = dsr[n] / d0;
            double ms0 = s[n] / (t + tw / 2 + tn / 2);
            Console.WriteLine(md0);
            Console.WriteLine(ms0);
            for (int k = 1; k < n; k++)
            {
                x[0] = 1;
                y[0] = 1;
                x[k] = x[k - 1] * ms[k];
                y[k] = y[k - 1] * md[k];
                h[k] = 0.25 * dsr[k] * (1 / y[k] * y[k] - 1 - 2.28 * rp[k] / dsr[k] + 0.56 * (rp[k] / dsr[k]) * (rp[k] / dsr[k])) / x[k] + rp[k];
                Console.WriteLine(h[k]);
            }
            double sd = t / d0;
            Console.WriteLine(sd);
        }


        // непонятно с циклами в массивах!
        public static void dpSmirAl()
        {
            Console.WriteLine("Введи диам. пуанс. свертки");
            double dpsw = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи диам. пуанс. окончат.выт ");
            double dpow = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи количество вытяжек ");
            int n = int.Parse(Console.ReadLine());

            double[] dp = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    dp[i] = dpsw;
                }
                else if (i == n)
                {
                    dp[i] = dpow;
                }
                else if (i > 0 && i < n)
                {
                    dp[i] = dp[i - 1] / Math.Pow((dpsw / dpow), (1 / n));
                    dp[i] = Math.Round(dp[i], 2);
                }
                Console.WriteLine(dp[i]);
            }
        }

        public static void одноконМатрПрижим(string[] args)
        {
            double pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            double alf = double.Parse(Console.ReadLine());
            double al = pi * alf / 180;
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double mum = 0.05;

            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(pi * alf / 180) + 1);
            dk = Math.Round(dk, 1);

            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);

            double rw = 3 * s;
            double rws = rw + s / 2;

            double md12 = d1 / dk;
            double md11 = dk / d0;
            double psi = 1 - dk / d0;
            double fik = pi * (90 - alf) / 180;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);

            double a = Math.Log(1 / md11) - psi + s / (2 * rws * Math.Sqrt(md11));
            double b = 1 + mum * fik;
            double c = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * (a * b / (1 - 0.2 * mum * b * c / md1) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);

            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Math.PI * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки без утонения=" + pbu + " кг");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Math.PI * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                //конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k); ;
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Math.PI * d1 * s * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmzk = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки =" + pu + " кг");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmzk = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки =" + pk + " кг");
            }
            double ustKonus = (Math.Sqrt((20 * sd) * (20 * sd) * (1 - Math.Sin(al)) + Math.Sin(al) * Math.Sin(al)) - 20 * sd) / Math.Sin(al);
            //проверка на склакообразованиие в конусе матрицы
            if (md1 < ustKonus)
            {
                Console.WriteLine("Однокон. матр.с плоск.и КОНИЧ. прижимом");
            }
            else
            {
                Console.WriteLine("Одноконус. матрица с прижимом");
            }
            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf + " град");
        }

        public static double sigmaz(double mum, double sigmr, double sigms2, double al, double ms1)
        {
            return ((1 + mum * (1 - sigmr / (1.15 * sigms2)) / Math.Sin(al) - mum * Math.Log(1 / ms1) / Math.Sin(al)) * Math.Log(1 / ms1) + sigmr / (1.15 * sigms2) + Math.Sin(al) / 2) * 1.15 * sigms2;
        }

        public static double sigmas2(double sigmb, double md1, double ms1, double psir)
        {
            return sigmb * Math.Pow((1 - md1 * ms1) / psir, psir / (1 - psir));
        }

        public static double sigmazk(double mum, double msk, double al, double sigms2k)
        {
            return 1.15 * sigms2k * ((1 + mum * (1 - Math.Log(1 / msk)) / Math.Sin(al)) * Math.Log(1 / msk) + Math.Sin(al) / 2);
        }

        static void одноконМатрБезПриж()
        {
            double pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи угол конусности матрицы");
            double alf = double.Parse(Console.ReadLine());
            double al = pi * alf / 180;
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double mum = 0.05;

            double dk = Math.Round(0.9 * d0, 1);
            double hk = Math.Round((dk - dm1) / 2 / Math.Tan(al), 1);
            double rw = Math.Round(0.05 * d0, 1);
            double rws = rw + s / 2;
            double dkk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(pi * alf / 180) + 1);

            double md12 = d1 / dkk;
            double md11 = dkk / d0;
            double psi = 1 - dkk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir))) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(al)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms={0} кг / мм^2", sigms);
                Console.WriteLine("sigmr={0} кг / мм^2", sigmr);
                Console.WriteLine("усилие свертки={0} кг", pbu);
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                //конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms={0} кг / мм^2", sigms2);
                Console.WriteLine("sigmr={0} кг / мм^2", sigmz);
                Console.WriteLine("усилие свертки={0} кг", pu);
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms={0} кг / мм^2", sigms2k);
                Console.WriteLine("sigmr={0} кг / мм^2", sigmzk);
                Console.WriteLine("усилие свертки={0} кг", pk);
            }
            Console.WriteLine("Одноконус. матрица без прижима");
            Console.WriteLine("Диаметр входн.кромки конуса={0}", dk);
            Console.WriteLine("Радиус входн.кромки конуса={0}", rw);
            Console.WriteLine("Высота конуса={0}", hk);
            Console.WriteLine("Угол конуса={0}", alf);
        }
    }
}
