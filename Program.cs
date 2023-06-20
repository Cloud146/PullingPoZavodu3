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
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, угол верхн. конуса, угол нижн. конуса
        //На выход: Диам. Верхн. конуса, Диам. Нижн. Конуса, Высота вехрн.,нижн. конусов, Углы верхн.,нижн. конусов, Радиус входн. верх. конуса, Радиус сопр. конусов, Усилие свертки
        public static void двухконМатр()
        {
            // Define constants
            double mum = 0.05;
            double Pi = Math.PI;

            // Read input values
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate recommended angle for upper cone
            if (sd > 0.012 && sd < 0.018)
            {
                string werhugol = "30 degrees";
                Console.WriteLine("Recommended angle for upper cone = " + werhugol);
            }
            else if (sd > 0.018 && sd < 0.05)
            {
                string werhugol = "45 degrees";
                Console.WriteLine("Recommended angle for upper cone = " + werhugol);
            }

            double alfw = double.Parse(Console.ReadLine());
            double alfn = double.Parse(Console.ReadLine());

            // Convert angles to radians
            double alw = Pi * alfw / 180;
            double alf = Pi * alfn / 180;

            double dk = d1 * Math.Sqrt((1 / (md1 * md1) - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double dw = 0.9 * d0;
            dw = Math.Round(dw, 1);
            double hw = (dw - dk) / (2 * Math.Tan(alw));
            hw = Math.Round(hw, 1);
            double rw = 0.05 * d0;
            rw = Math.Round(rw, 1);
            double rs = (d0 - dm1) / 3;
            rs = Math.Round(rs, 1);
            double r = (d0 - d1) / 2;
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double fik = (alw - alf) / 2;
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (((1 + mum / Math.Tan(alw)) * (Math.Log(1 / md11) - psi) + s / (2 * d1 * Math.Sqrt(md11))) * (1 + mum * fik) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);

            if (ms1 == 1)
            {
                // усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("Усилие свертки = " + pbu + " kg"); // Только для консоли
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("Final stage");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("Force of rolling = " + pk + " kg");
            }
            Console.WriteLine("sigms = " + sigms + " kg / mm^2");
            Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");

            Console.WriteLine("Двухконус. матрица");
            Console.WriteLine("Диам. верхнего конуса", dw);
            Console.WriteLine("Диаметр нижнего конуса", dk);
            Console.WriteLine("Высота верхнего конуса", hw);
            Console.WriteLine("Высота нижнего конуса", hk);
            Console.WriteLine("Угол верхнего конуса", alfw);
            Console.WriteLine("угол нижнего конуса", alfn);
            Console.WriteLine("Радиус входн.верх конуса", rw);
            Console.WriteLine("Радиус сопряж. конусов", rs);
        }

        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ
        //На выход: sigmas, sigmaz, Усилие свертки, Усилие свертки в конечной стадии, Радиус матрицы, Тип прижима(Плоский или плоский и торцевой)
        public static void радМатрПриж()
        {
            // Define the constants.
            const double mum = 0.05;
            const double Pi = Math.PI;

            // Get the values from the Excel spreadsheet.
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate the radius of the matrix.
            double rms = 0.69 * d1 * (0.16 * Math.Sqrt(18.3 / md1 * md1 + 21.2 + 10.2 * (rps / d1) * (rps / d1) - 41.3 * rps / d1) - 1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);

            // Calculate the reduced modulus of deformation.
            double md11 = (d1 + 2 * rms) / d0;
            double md12 = md1 / md11;
            double alr = Math.Atan2(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)), rm + s1);
            double fi = Pi / 2 - alr;
            double al = alr / 2;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * (psisr / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);

            // Calculate the shear modulus.
            double a = Math.Log(1 / md1) + s / (4 * rm);
            double b = 1 + 1.5 * mum;
            double c = 0.2 * mum * b / md1;
            double d = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * a * b / (1 - c * d);
            sigmr = Math.Round(sigmr, 1);

            // Calculate the rolling force without thinning.
            double pbu;

            // Calculate the rolling force with thinning.
            if (ms1 == 1)
            {
                //усилие при свертке без утонения.
                pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                //конечная стадия
                // Calculate the rolling force with thinning in the final stage.
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }

            Console.WriteLine("Радиальная матрица с прижимом");
            if (sd < (1 - psir * Math.Log(1 / md1) - md1) / 20)
            {
                Console.WriteLine("Требуется плоский и торцевой прижим");
            }
            else
            {
                Console.WriteLine("Плоский прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }

        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ
        //На выход: sigmas, sigmaz, Усилие свертки, Усилие свертки в конечной стадии, Радиус матрицы
        public static void радМатрБезПриж()
        {
            // Define the variables.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Calculate the radius of the matrix.
            double rms = d1 * (1 - md1) / (2 * md1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);

            // Calculate the other parameters.
            double fig = 200 * (1 - md1 * (1 + 2.28 * rps / d1)) / (1 + md1 * (5 - 6 * md1));
            double fi = Pi * fig / 180;
            double md11 = md1 * (1 + 2 * rms * (1 - Math.Cos(fi)) / d1);
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double al = alr / 2;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * (psisr / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum * fi) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);

            // Calculate the required forces.
            double pbu, pu;
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }

            // Print the results.
            Console.WriteLine("Радиальная матрица с прижимом");
            if (sd < (1 - psir * Math.Log(1 / md1) - md1) / 20)
            {
                Console.WriteLine("Требуется плоский и торцевой прижим");
            }
            else
            {
                Console.WriteLine("Плоский прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }
        //На вход: Тип матрицы (Коническая, Радиальная), Тип матрицы (Однокон., Двухкон., Однокон. с радиус.), sd, md1, psir
        //На выход: РАЗВЕТВЛЕНИЕ ПО МЕТОДАМ
        public static void ВырСвертка(string[] args)
        {

            // Get the values from the user
            Console.WriteLine("Введите тип матрицы: матрица конич.-K, радиальная-R");
            string tipMatriz = Console.ReadLine();
            double sd = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Check if the squeeze is needed
            if (sd > (1 - md1) / 18)
            {
                Console.WriteLine("прижим не нужен на кон. и рад. матр");
            }
            else if (sd > (1 - md1) / 36)
            {
                Console.WriteLine("на конич. матр. прижим не нужен");
            }
            else
            {
                Console.WriteLine("нужен прижим");
            }

            // Check the type of matrix
            if (tipMatriz == "k")
            {
                if (sd > 0.012 && sd < 0.05)
                {
                    Console.WriteLine("Рекомендуется двухконусная матрица");
                }
                else if (sd > 0.05)
                {
                    Console.WriteLine("Рекомендуется одноконусная матрица с большим радиусом");
                }

                // Check if the squeeze is needed for conical matrix
                if (sd <= (1 - md1) / 36)
                {
                    Console.WriteLine("одноконМатрПрижим");
                }
                else
                {
                    Console.WriteLine("матрица однокон.-O, двухкон.-D, однокон. с радиус-OR");
                    string tipKonMatriz = Console.ReadLine();

                    if (tipKonMatriz == "o")
                    {
                        Console.WriteLine("одноконМатрБезПриж");
                    }
                    else if (tipKonMatriz == "d")
                    {
                        Console.WriteLine("двухконМатр");
                    }
                    else if (tipKonMatriz == "or")
                    {
                        Console.WriteLine("одноконМатрРадиус");
                    }
                }
            }
            else if (tipMatriz == "r")
            {
                if (sd > (1 - md1) / 18)
                {
                    Console.WriteLine("радМатрБезПриж");
                }
                else
                {
                    Console.WriteLine("радМатрПрижим");
                }
            }
        }

        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ
        //На выход: sigmas, sigmaz, Усилие свертки, Усилие свертки в конечной стадии, Диам. Конуса, Высота конуса, Угол Конуса, Радиус Конуса
        public static void одноконМатрРадиус(string[] args)
        {
            // Define constants
            double mum = 0.05;
            double Pi = Math.PI;

            // Get input values
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            //Введи угол  конуса
            double alfn = double.Parse(Console.ReadLine());
            double alf = Pi * alfn / 180;

            // Calculate values
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double md11 = dk / d0;
            double a = md11 + sd;
            double b = (1 - Math.Sin(alf)) * Math.Tan(alf);
            double hm = (a * (1 - b) + b - sd * ms1 - md1) * d0 / (2 * Math.Tan(alf));
            hm = Math.Round(hm, 1);
            double c = hm - hk;
            if (c < 0)
            {
                Console.WriteLine("высота матр. меньше высоты конуса");
                return;
            }
            double dkm = dm1 + 2 * hm * Math.Tan(alf);
            dkm = Math.Round(dkm, 1);
            double md12 = d1 / dk;
            double rw = (d0 - dk - s) / 2;
            rw = Math.Round(rw, 1);
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir))) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(alf)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            if (ms1 == 1)
            {
                // усилие при свертке без утонения
                double pbu = Math.PI * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                //усилие при свертке с утонением
                double pu = Pi * d1 * s1 * sigmz;
                pu = Math.Round(pu, 0);
                // конечная стадия
                double msk = ms1 * Math.Sqrt(md1);
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);

                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("Диам.  конуса=" + dkm);
            Console.WriteLine("Высота конуса=" + hm);
            Console.WriteLine("угол  конуса=" + alfn);
            Console.WriteLine("радиус конуса=" + rw);
        }

        //На вход: Количество матриц
        //На выход: Углы конустности матриц
        public static void УголКонусМатр(string[] args)
        {
            // Define constants
            double Pi = Math.PI;

            // Get input values
            int n = int.Parse(Console.ReadLine());
            double mum = 0.05;

            // Create arrays to store the input values
            double[] s = new double[n - 1];
            double[] d = new double[n - 1];
            double[] md = new double[n - 1];

            // Get the input values and store them in the arrays
            for (int i = 0; i < n - 1; i++)
            {
                s[i] = double.Parse(Console.ReadLine());
                d[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the angles
            for (int i = 0; i < n - 1; i++)
            {
                double a = Math.Sqrt(0.76 * mum * Math.Log(1 / md[i]) * Math.Sqrt(d[i] / s[i]));
                double b = Math.Asin(a);
                double alf = 180 * b / Pi;
                alf = Math.Round(alf, 0);
                Console.WriteLine("alf={0}", alf);
            }
        }
        //На вход: Количество матриц
        //На выход: Высоты конустности матриц
        public static void ВысотаКонусаМатр(string[] args)
        {
            // Define constants
            double Pi = Math.PI;

            // Get input values
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the input values
            double[] dm = new double[n - 1];
            double[] alf = new double[n - 1];

            // Get the input values and store them in the arrays
            for (int i = 0; i < n - 1; i++)
            {
                dm[i] = double.Parse(Console.ReadLine());
            }

            // Get the input radius
            //Введи радиус вход. конуса матр
            double r = double.Parse(Console.ReadLine());

            // Calculate the heights
            for (int i = 1; i < n - 1; i++)
            {
                alf[i] = double.Parse(Console.ReadLine());
                double al = Pi * alf[i] / 180;
                double hm = (dm[i - 1] - dm[i]) / 2 / Math.Tan(al) + r * (1 - Math.Sin(al));
                double dkm = dm[i - 1] + 2 * r * (1 - Math.Sin(al)) * Math.Tan(al);
                dkm = Math.Round(dkm, 1);
                hm = Math.Round(hm, 1);
                Console.WriteLine("hm={0}", hm);
                Console.WriteLine("dkm={0}", dkm);
                Console.WriteLine("r={0}", r);
            }
        }
        //На вход: md, sigmamb, psir
        //На выход: sigms
        public static void ВытБезУтон()
        {
            // Define constants
            double Pi = Math.PI;
            double md = double.Parse(Console.ReadLine());
            double sigmb = double.Parse(Console.ReadLine());
            double psisr = double.Parse(Console.ReadLine());
            // Calculate the stress
            double psir = 1 - Math.Sqrt(md);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir))) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            Console.WriteLine("sigms={0}", sigms);
        }
        //На вход: Массивы alf, md, ms, s, dsr
        //На выход: sigmz, sigms2, p и sigms2k, sigmzk, pk
        static void УсилиеВытяжки()
        {
            // Define constants
            const double Pi = Math.PI;
            const double mum = 0.05;

            // Get the number of iterations
            int n = int.Parse(Console.ReadLine());

            // Get the material properties
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Create arrays to store the data
            double[] alf = new double[n - 1];
            double[] md = new double[n - 1];
            double[] ms = new double[n - 1];
            double[] s = new double[n];
            double[] dsr = new double[n];

            // Read the data from the console
            for (int i = 0; i < n - 1; i++)
            {
                alf[i] = double.Parse(Console.ReadLine());
                md[i] = double.Parse(Console.ReadLine());
                ms[i] = double.Parse(Console.ReadLine());
            }
            for (int i = 0; i < n; i++)
            {
                s[i] = double.Parse(Console.ReadLine());
                dsr[i] = double.Parse(Console.ReadLine());
            }

            // Calculate the results
            for (int i = 0; i < n - 1; i++)
            {
                // Convert the angle from degrees to radians
                double alfaRad = alf[i] * Pi / 180;

                if (ms[i] < 1)
                {
                    // Calculate the yield stresses
                    double sigms1 = sigmas1(sigmb, md[i], psir);
                    double sigms2 = sigmas2(sigmb, md[i], ms[i], psir);
                    sigms2 = Math.Round(sigms2, 1);

                    // Calculate the reduced yield stresses
                    double sigmr1 = 1.1 * sigms1 * sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);
                    double sigr1 = sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);

                    // Calculate the bursting stresses
                    double sigmz = sigmaz2(ms[i], alfaRad, sigr1, sigms2);
                    sigmz = Math.Round(sigmz, 1);

                    // Calculate the bursting forces
                    double pst = pste(dsr[i], s[i], sigmz);
                    double ptr = Pi * ptre(ms[i], dsr[i], s[i], alfaRad, sigms2, sigr1);
                    double p = pst + ptr;
                    p = Math.Round(p, 0);

                    //конец вытяжки с утонением
                    double msk = s[i] * Math.Sqrt(md[i]) / s[i - 1];
                    double sigms2k = sigmas2(sigmb, md[i], msk, psir);
                    sigms2k = Math.Round(sigms2k, 1);
                    double sigmr1k = 0;
                    double sigmzk = sigmaz2(msk, alfaRad, sigr1, sigms2k);
                    sigmzk = Math.Round(sigmzk, 1);
                    double pstk = pste(dsr[i], s[i], sigmzk);
                    double ptrk = Pi * ptre(msk, dsr[i], s[i], alfaRad, sigms2k, sigr1);
                    double pk = pstk + ptrk;
                    pk = Math.Round(pk, 0);
                    Console.WriteLine("p = " + p);
                    Console.WriteLine("sigmz = " + sigmz);
                    Console.WriteLine("sigms2 = " + sigms2);
                    Console.WriteLine("pk = " + pk);
                    Console.WriteLine("sigmzk = " + sigmzk);
                    Console.WriteLine("sigms2k = " + sigms2k);
                }
                else if (ms[i] == 1)
                {
                    double psisr = 1 - Math.Pow(md[i], 2);
                    double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
                    double sigmr = 1.1 * sigms * sigmar1(alfaRad, md[i], s[i - 1], dsr[i - 1]);
                    double p = Pi * dsr[i] * s[i] * sigmr;
                    sigms = Math.Round(sigms, 1);
                    sigmr = Math.Round(sigmr, 1);
                    p = Math.Round(p, 0);
                    Console.WriteLine("p = " + p);
                    Console.WriteLine("sigmr = " + sigmr);
                    Console.WriteLine("sigms = " + sigms);
                }
            }
        }
        //На вход: alf, md, s, dsr
        //На выход: sigmar1
        public static double sigmar1(double alf, double md, double s, double dsr)
        {
            double mum = 0.05;
            return (1 + mum / Math.Tan(alf)) * Math.Log(1 / md) + 0.66 * Math.Sin(alf) * Math.Sqrt(s / dsr);
        }
        //На вход: double ms, double alf, double sigmr1, double sigms2
        //На выход: sigmaz2
        public static double sigmaz2(double ms, double alf, double sigmr1, double sigms2)
        {
            double mum = 0.05;
            double mup = 0.05;
            double ks = 1 / ms;
            double a = 1 - 0.5 * Math.Log(ks) - sigmr1;
            double b = 0.5 * (Math.Log(ks)) * (Math.Log(ks)) - sigmr1 * (ks - 1 - Math.Log(ks));
            return 1.15 * sigms2 * (Math.Log(ks) + mum * a * Math.Log(ks) / Math.Sin(alf) - mup * b / Math.Sin(alf) + sigmr1 + (1 - Math.Cos(alf)) / Math.Sin(alf));
        }
        //На вход: d, s, sigmz
        //На выход: pste
        public static double pste(double d, double s, double sigmz)
        {
            double Pi = Math.PI;
            return Pi * d * s * sigmz;
        }
        //На вход: double ms, double d, double s, double alf, double sigms2, double sigmr1
        //На выход: ptre
        public static double ptre(double ms, double d, double s, double alf, double sigms2, double sigmr1)
        {
            double ks = 1 / ms;
            double a = Math.Log(ks) - sigmr1 * (ks - 1);
            return 1.15 * 0.05 * sigms2 * d * s * a / Math.Sin(alf);
        }
        //На вход: sigmb, md, psir
        //На выход: sigmas1
        public static double sigmas1(double sigmb, double md, double psir)
        {
            return sigmb * Math.Pow(((1 - md) / psir), (psir / (1 - psir)));
        }
        //На вход:  narDiam,  dopuNarZiam,  vnutrDiam, dopuVnutrDiam, vysota,  tolshchLenty,  verhDopTolst,  nizDopTolst,  radius
        //На выход: D0KolpShtBezUt
        public static double D0KolpShtBezUt(double narDiam, double dopuNarZiam, double vnutrDiam, double dopuVnutrDiam, double vysota, double tolshchLenty, double verhDopTolst, double nizDopTolst, double radius)
        {
            double dn = narDiam + dopuNarZiam / 2;
            double dw = vnutrDiam + dopuVnutrDiam / 2;
            double h = vysota;
            double s = tolshchLenty + (verhDopTolst + nizDopTolst) / 2;
            double r = radius;
            double Pi = Math.PI;
            double vcil = Pi * (dn + dw) * s * (h - s - r) / 2;
            double vdno = Pi * s * (dw - 2 * r) * (dw - 2 * r) / 4;
            double vtor = Pi * Pi * (dw / 2 - r) * ((r + s) * (r + s) - r * r) / 2;
            double v = vcil + vdno + vtor;
            double d0KolpShtBezUt = Math.Sqrt(4 * v / Pi / s);
            return Math.Round(d0KolpShtBezUt, 2);
        }
        //На вход: md1, sd, psir, тип матрицы(Коническая или Радиальная) и Двухкон. или Однокон.
        //На выход: РАЗВЕТВЛЕНИЕ и ТИП ПРИЖИМА
        static void вырСверт2матр()
        {
            // Get the values from the user.
            Console.WriteLine("Enter the value of md1:");
            double md1 = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter the value of sd:");
            double sd = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter the value of psir:");
            double psir = Convert.ToDouble(Console.ReadLine());

            // Check if the shear stress is greater than (1 - md1) / 18.
            if (sd > (1 - md1) / 18)
            {
                Console.WriteLine("прижим не нужен на кон. и рад. матр");
            }
            else if (sd > (1 - md1) / 36)
            {
                Console.WriteLine("на конич. матр. прижим не нужен");
            }
            else
            {
                Console.WriteLine("нужен прижим");
            }

            // Get the type of matrix from the user.
            Console.WriteLine("Enter the type of matrix (матрица конич.-K, радиальная-R):");
            string tipMatriz = Console.ReadLine();

            // Check if the matrix type is conical.
            if (tipMatriz == "k")
            {
                if (sd > 0.012 && sd < 0.05)
                {
                    Console.WriteLine("Рекомендуется двухконусная матрица");
                }
                else if (sd > 0.05)
                {
                    Console.WriteLine("Рекомендуется одноконусная матрица с большим радиусом");
                }

                // Check if the shear stress is less than or equal to (1 - md1) / 36.
                if (sd <= (1 - md1) / 36)
                {
                    одноконМатрПрижим2матр();
                }
                else
                {
                    Console.WriteLine("матрица однокон.-O, двухкон.-D, однокон. с радиус-OR");
                    string tipKonMatriz = Console.ReadLine();
                    if (tipKonMatriz == "o")
                    {
                        одноконМатрБезПриж2матр();
                    }
                    else if (tipKonMatriz == "d")
                    {
                        двухконМатр2матр();
                    }
                    else if (tipKonMatriz == "or")
                    {
                        одноконМатрРадиус2матр();
                    }
                }
            }
            // Check if the matrix type is radial.
            else if (tipMatriz == "r")
            {
                if (sd > (1 - md1) / 18)
                {
                    радМатрБезПриж2матр();
                }
                else
                {
                    радМатрПрижим2матр();
                }
            }
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, Угол Конуст. Верхней матрицы, Коеф. утонения на верх. матр.
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        // Диаметр входн.кромки конуса, Радиус входн.кромки конуса, Высота конуса, Угол конуса
        public static void одноконМатрПрижим2матр()
        {
            // Get the input values.
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Get the user input for the angle of cone of the upper matrix.
            Console.WriteLine("Введи угол конусности верх. матрицы");
            double alf = double.Parse(Console.ReadLine());
            double al = Pi * alf / 180;

            double ms1 = double.Parse(Console.ReadLine());
            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the coefficient of friction.
            double mum = 0.05;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(Math.PI * alf / 180) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);
            double rw = 3 * s;
            double rws = rw + s / 2;

            // Calculate the force on the upper matrix
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double psi = 1 - dk / d0;
            double fik = Math.PI * (90 - alf) / 180;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double a = Math.Log(1 / md11) - psi + s / (2 * rws * Math.Pow(md11, 2));
            double b = 1 + mum * fik;
            double c = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * (a * b / (1 - 0.2 * mum * b * c / md1) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * (psisrn / psir) * Math.Pow(psir / (1 - psir), 1 / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1sr == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s1 * sigmzk;
                //усилие при свертке с утонением в конечн. стадии;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Одноконус.с приж.через 2 матр.");
            double ustKonus = (Math.Sqrt((20 * sd) * (20 * sd) * (1 - Math.Sin(al)) + Math.Sin(al) * Math.Sin(al)) - 20 * sd) / Math.Sin(al);
            //проверка на склакообразованиие в конусе матрицы
            if (md1 < ustKonus)
            {
                Console.WriteLine("Однокон.с плоск.и КОНИЧ. приж. через 2 матр.");
            }

            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf);
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, Угол Конуст. Верхней матрицы, Коеф. утонения на верх. матр.
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        // Диаметр входн.кромки конуса, Радиус входн.кромки конуса, Высота конуса, Угол конуса
        public static void одноконМатрБезПриж2матр()
        {
            // Get the input values.
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            // Get the user input for the angle of cone of the upper matrix.
            Console.WriteLine("Введи угол конусности верх. матрицы");
            double alf = double.Parse(Console.ReadLine());
            double al = Pi * alf / 180;

            double ms1 = double.Parse(Console.ReadLine());
            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the coefficient of friction.
            double mum = 0.05;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = 0.9 * d0;
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(al);
            hk = Math.Round(hk, 1);
            double rw = 0.5 * d0;
            double rws = rw + s / 2;
            double dkk = d1 * Math.Sqrt((1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.07 * alf * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(Math.PI * alf / 180) + 1);


            // Calculate the force on the upper matrix
            double md12 = d1 / dkk;
            double md11 = dkk / d0;
            double psi = 1 - dkk / d0;
            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);

            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(al)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                double pk = Pi * d1 * s1 * sigmzk;
                //усилие при свертке с утонением в конечн. стадии;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Одноконус.матр.без приж.через 2 матр");
            Console.WriteLine("Диаметр входн.кромки конуса = " + dk);
            Console.WriteLine("Радиус входн.кромки конуса = " + rw);
            Console.WriteLine("Высота конуса = " + hk);
            Console.WriteLine("Угол конуса = " + alf);
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, Средний коеф. утонения, Коеф. утонения на верх. матр., 
        //радиус вход. конуса нижн. матр, угол верх. конуса, угол нижнего конуса
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        // Диаметры конусов, Радиус входн.кромки конуса, Высоты конусов, Углы конусов, Радиус сопряж. конусов
        public static void двухконМатр2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());

            // Calculate the average coefficient of thinning.
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            // Get the user input for the coefficient of thinning of the upper matrix.
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            // Calculate the coefficient of thinning of the lower matrix.
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfanm = Math.Asin(q);
            double alfnm = alfanm * 180 / Pi;
            alfnm = Math.Round(alfnm, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfanm) + rn * (1 - Math.Sin(alfanm));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfanm)) * Math.Tan(alfanm);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //геометрия верхней матрицы
            double sd = double.Parse(Console.ReadLine());
            double werhugol;
            if (sd > 0.012 && sd < 0.018)
            {
                werhugol = 30;
                Console.WriteLine("Рекомендуемый угол верх. конуса = " + werhugol);
            }
            else if (sd > 0.0018 && sd < 0.05)
            {
                werhugol = 45;
                Console.WriteLine("Рекомендуемый угол верх. конуса = " + werhugol);
            }
            Console.WriteLine("Введи угол верх. конуса");
            double alfw = double.Parse(Console.ReadLine());
            Console.WriteLine("Введи угол нижнего. конуса");
            double alfn = double.Parse(Console.ReadLine());
            double alw = Pi * alfw / 180;
            double alf = Pi * alfn / 180;

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * (rps / d1) * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfanm) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfn * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double dw = 0.9 * d0;
            dw = Math.Round(dw, 1);
            double hw = (dw - dk) / (2 * Math.Tan(alw));
            hw = Math.Round(hw, 1);
            double rw = 0.5 * d0;
            rw = Math.Round(rw, 1);
            double rs = (d0 - dm1) / 3;
            rs = Math.Round(rs, 1);
            double r = (d0 - d1) / 5;

            // Calculate the force on the upper matrix
            double md12 = d1 / dk;
            double md11 = dk / d0;
            double fik = (alw + alf) / 2;
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (((1 + mum / Math.Tan(alw)) * (Math.Log(1 / md11) - psi) + s / (2 * r * Math.Sqrt(md11))) * (1 + mum * fik) + Math.Log(1 / md12));
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfanm / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfanm);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfnm);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Двухкон.матр2 матр");
            Console.WriteLine("Диам. верхнего конуса = " + dw);
            Console.WriteLine("Диаметр нижнего конуса = " + dk);
            Console.WriteLine("Высота верхнего конуса = " + hw);
            Console.WriteLine("Высота нижнего конуса = " + hk);
            Console.WriteLine("Угол верхнего конуса = " + alfw);
            Console.WriteLine("Угол верхнего конуса = " + alfn);
            Console.WriteLine("Радиус входн.верх конуса = " + rw);
            Console.WriteLine("Радиус входн.верх конуса = " + rs);
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, Средний коеф. утонения, Коеф. утонения на верх. матр., 
        //радиус вход. конуса нижн. матр
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        // Диаметр входн.кромки конуса, Радиус входн.кромки конуса, Высота конуса, Угол конуса
        public static void одноконМатрРадиус2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double s1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);

            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());

            double ms1n;
            //нижняя матрица
            if (ms1w == 1)
            {
                ms1n = ms1;
            }
            else
            {
                ms1n = ms1 / ms1w;
                ms1n = Math.Round(ms1n, 3);
            }

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            s1 = s * ms1;

            // Get the user input for the diameter of the punch.
            double dp = double.Parse(Console.ReadLine());

            // Calculate the diameter of the upper matrix.
            double dm1;
            if (ms1w == 1)
            {
                dm1 = dp + 1.5 * s1 / Math.Sqrt(md1);
                d1 = (dm1 + dp) / 2;
            }
            else
            {
                d1 = dp + s1;
                dm1 = dp + 2 * s1;
            }
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;

            // Calculate the angle of cone of the lower matrix.
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            // Get the user input for the radius of the inlet cone of the lower matrix.
            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            Console.WriteLine("Введи угол конуса верх. матрицы");
            double alfa = double.Parse(Console.ReadLine());
            double alf = Pi * alfa / 180;

            //расстояние между матрицами
            // Calculate the distance between the matrices.
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            // Calculate the geometry of the upper matrix
            double dk = d1 * Math.Sqrt((1 / md1 * md1 - 1 - 2.28 * rps / d1 + 0.07 * alfa * rps / d1 + 0.56 * (rps / d1) * (rps / d1)) * Math.Sin(alf) + 1);
            dk = Math.Round(dk, 1);
            double hk = (dk - dm1) / 2 / Math.Tan(alf);
            hk = Math.Round(hk, 1);
            double md11 = dk / d0;
            double a = md11 + sd;
            double b = (1 - Math.Sin(alf)) * Math.Tan(alf);
            double hm = (a * (1 - b) + b - sd * ms1 - md1) * d0 / (2 * Math.Tan(alf));
            hm = Math.Round(hm, 1);
            double c = hm - hk;
            if (c < 0)
            {
                Console.WriteLine("высота матр. меньше высоты конуса");
                return;
            }
            double dkm = dm1 + 2 * hm * Math.Tan(alf);
            dkm = Math.Round(dkm, 1);
            double md12 = d1 / dk;
            double rw = (d0 - dk - s) / 2;
            rw = Math.Round(rw, 1);
            double psi = 1 - dk / d0;
            double psisr = 1 - Math.Sqrt(md1);

            double sigms = sigmb * Math.Pow((psisr / psir), (psir / (1 - psir)) / (1 - psir));
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum / Math.Tan(alf)) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            //усилие вытяжки на нижней матрице
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow((psisrn / psir), (psir / (1 - psir)) / (1 - psir));
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Pi * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);
            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Однокон.матр.с больш.радиусом через 2матр");
            Console.WriteLine("Диам. конуса = " + dkm);
            Console.WriteLine("Высота конуса = " + hm);
            Console.WriteLine("Угол  конуса = " + alfa);
            Console.WriteLine("Радиус конуса = " + rw);
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, радиус вход. конуса нижн. матр, коэф. утонен. на верх. матр
        //радиус вход. конуса нижн. матр, Средний коэф. утонен
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр, Радиус МАТРИЦЫ
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        public static void радМатрПрижим2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());
            //нижняя матрица
            double ms1n = ms1 / ms1w;
            ms1n = Math.Round(ms1n, 3);

            // Set the coefficient of thinning of the upper matrix.
            ms1 = ms1w;
            double s1 = s * ms1;
            double dp = double.Parse(Console.ReadLine());
            d1 = dp + s1;
            //диаметр верхней матрицы
            double dm1 = dp + 2 * s1;
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;
            // угол кон. ниж. матр
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());

            // Calculate the height of the inlet cone of the lower matrix.
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //верхней матрицы
            double rms = 0.69 * d1 * (0.16 * Math.Sqrt(18.3 / md1 * md1 + 21.2 + 10.2 * (rps / d1) * (rps / d1) - 41.3 * rps / d1) - 1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);
            double md11 = (d1 + 2 * rms) / d0;
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double fi = Math.PI / 2 - alr;
            double al = alr / 2;

            //расстояние между матрицами
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret;
            if (ms1w == 1)
            {
                tret = 0;
            }
            else
            {
                tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(al));
            }
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double a = Math.Log(1 / md1) + s / (4 * rm);
            double b = 1 + 1.5 * mum;
            double c = 0.2 * mum * b / md1;
            double d = 1 - 18 * sd / (1 - md1);
            double sigmr = 1.1 * sigms * a * b / (1 - c * d);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow(psisrn / psir, psir / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);

            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, al, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, al, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Радиал. матрица с прижимом 2 матр");
            if (sd < ((md11 - md1) / 20))
            {
                Console.WriteLine("Треб.плоск.и тор.прижим");
            }
            else
            {
                Console.WriteLine("Треб.плоск.и тор.прижим");
            }
            Console.WriteLine("Радиус матрицы = " + rm);
        }
        //На вход: нарДиам, допуНарДиам, внутрДиам, допуВнутрДиам, высота, радиус, толщЛенты, верхДопТолщ, нижДопТолщ, радиус вход. конуса нижн. матр, коэф. утонен. на верх. матр
        //радиус вход. конуса нижн. матр, Средний коэф. утонен
        //На выход: усилие свертки, конечная стадия, sigms, sigmz, диам.верх. матр, коэф. утон.верх. матр, коэф. утон.ниж. матр, Радиус МАТРИЦЫ
        //угол конус.ниж.матр, радиус конус.ниж.матр, диам. конус.ниж.матр, высота конус.ниж.матр, расст. между поясками матр, sigmsn, sigmzn, усилие на нижней матр
        public static void радМатрБезПриж2матр()
        {
            // Get the input values.
            double mum = 0.05;
            double Pi = Math.PI;
            double d0 = double.Parse(Console.ReadLine());
            double d1 = double.Parse(Console.ReadLine());
            double md1 = double.Parse(Console.ReadLine());
            double dm1n = double.Parse(Console.ReadLine());
            double rp = double.Parse(Console.ReadLine());
            double s = double.Parse(Console.ReadLine());
            double ms1 = double.Parse(Console.ReadLine());
            double ms1sr = Math.Sqrt(ms1);
            Console.WriteLine("Средний коэф. утонен=" + ms1sr);
            Console.WriteLine("Введи коэф. утонен. на верх. матр");
            double ms1w = double.Parse(Console.ReadLine());
            double rps = rp + s / 2;
            double sigmb = double.Parse(Console.ReadLine());
            double psir = double.Parse(Console.ReadLine());
            double sd = double.Parse(Console.ReadLine());

            //нижняя матрица
            double ms1n = ms1 / ms1w;
            ms1n = Math.Round(ms1n, 3);
            ms1 = ms1w;
            double s1 = s * ms1;
            double dp = double.Parse(Console.ReadLine());
            d1 = dp + s1;
            //диаметр верхней матрицы
            double dm1 = dp + 2 * s1;
            dm1 = Math.Round(dm1, 2);
            md1 = d1 / d0;
            // угол кон. ниж. матр
            double q = Math.Sqrt(2 * mum * Math.Log(1 / ms1n) * (1 - Math.Log(1 / ms1n)));
            double alfan = Math.Asin(q);
            double alfn = alfan * 180 / Pi;
            alfn = Math.Round(alfn, 0);

            Console.WriteLine("Введи радиус вход. конуса нижн. матр");
            double rn = double.Parse(Console.ReadLine());
            double hmn = (dm1 - dm1n) / 2 / Math.Tan(alfan) + rn * (1 - Math.Sin(alfan));
            double dkmn = dm1 + 2 * rn * (1 - Math.Sin(alfan)) * Math.Tan(alfan);
            dkmn = Math.Round(dkmn, 1);
            hmn = Math.Round(hmn, 1);

            //верхней матрицы
            double rms = d1 * (1 - md1) / (2 * md1);
            double rm = rms - s1 / 2;
            rm = Math.Round(rm, 1);
            double fig = 200 * (1 - md1 * (1 + 2.28 * rps / d1)) / (1 + md1 * (5 - 6 * md1));
            double fi = fig * Pi / 180;
            double md11 = md1 * (1 + 2 * rms * (1 - Math.Cos(fi)) / d1);
            double md12 = md1 / md11;
            double alr = Math.Atan(Math.Sqrt((rm + s) * (rm + s) - (rm + s1) * (rm + s1)) / (rm + s1));
            double alf = alr / 2;

            //расстояние между матрицами
            double perw = 1 / Math.Pow(md1, 2) - 1 - 2.28 * rps / d1 + 0.56 * Math.Pow(rps / d1, 2);
            double wtor = Math.Sqrt(1 - Math.Pow(ms1w, 2) + 2 * rps * (1 - ms1w) / s);
            double tret;
            if (ms1w == 1)
            {
                tret = 0;
            }
            else
            {
                tret = ms1w * (1 - ms1n) / Math.Tan(alfan) - (1 / md1 - Math.Pow(ms1w, 2)) / (4 * ms1w * Math.Tan(alf));
            }
            double rast = 0.25 * d1 * perw / ms1w + s * (wtor + tret);
            rast = Math.Round(rast, 1);

            double psisr = 1 - Math.Sqrt(md1);
            double sigms = sigmb * Math.Pow(psisr / psir, psir / (1 - psir)) / (1 - psir);
            sigms = Math.Round(sigms, 1);
            double sigmr = 1.1 * sigms * (1 + mum * fig) * Math.Log(1 / md12);
            sigmr = Math.Round(sigmr, 1);
            double psisrn = 1 - 0.5 * md1 * ms1w * (1 + ms1n);
            double sigmsn = sigmb * Math.Pow(psisrn / psir, psir / (1 - psir)) / (1 - psir);
            double sigmzn = sigmsn * (Math.Log(1 / ms1n) + alfan / 2 + mum * (1 - Math.Log(1 / ms1n) * Math.Log(1 / ms1n)) / alfan);
            double sn = s * ms1w * ms1n;
            double pn = Math.PI * (dp + sn) * sn * sigmzn;
            pn = Math.Round(pn, 0);
            sigmsn = Math.Round(sigmsn, 1);
            sigmzn = Math.Round(sigmzn, 1);

            if (ms1 == 1)
            {
                //усилие при свертке без утонения
                double pbu = Pi * d1 * s * sigmr;
                pbu = Math.Round(pbu, 0);
                Console.WriteLine("sigms = " + sigms + " kg / mm^2");
                Console.WriteLine("sigmr = " + sigmr + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pbu + " kg");
            }
            else
            {
                double sigms2 = sigmas2(sigmb, md1, ms1, psir);
                sigms2 = Math.Round(sigms2, 1);
                double sigmz = sigmaz(mum, sigmr, sigms2, alf, ms1);
                sigmz = Math.Round(sigmz, 1);
                double pu = Pi * d1 * s1 * sigmz;
                //усилие при свертке с утонением;
                pu = Math.Round(pu, 0);
                double msk = ms1 * Math.Sqrt(md1);
                //конечная стадия;
                double sigmrk = 0;
                double sigms2k = sigmas2(sigmb, md1, msk, psir);
                sigms2k = Math.Round(sigms2k, 1);
                double sigmzk = sigmazk(mum, msk, alf, sigms2k);
                sigmzk = Math.Round(sigmzk, 1);
                //усилие при свертке с утонением в конечн. стадии;
                double pk = Pi * d1 * s1 * sigmzk;
                pk = Math.Round(pk, 0);
                Console.WriteLine("sigms = " + sigms2 + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmz + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pu + " kg");
                Console.WriteLine("конечная стадия");
                Console.WriteLine("sigms = " + sigms2k + " kg / mm^2");
                Console.WriteLine("sigmz = " + sigmzk + " kg / mm^2");
                Console.WriteLine("усилие свертки = " + pk + " kg");
            }
            Console.WriteLine("диам.верх. матр.= " + dm1);
            Console.WriteLine("коэф. утон.верх. матр = " + ms1w);
            Console.WriteLine("коэф. утон.ниж. матр = " + ms1n);
            Console.WriteLine("угол конус.ниж.матр = " + alfn);
            Console.WriteLine("радиус конус.ниж.матр = " + rn);
            Console.WriteLine("диам. конус.ниж.матр = " + dkmn);
            Console.WriteLine("высота конус.ниж.матр = " + hmn);
            Console.WriteLine("расст. между поясками матр = " + rast);
            Console.WriteLine("sigmsn = " + sigmsn + " kg / mm^2");
            Console.WriteLine("sigmzn = " + sigmzn + " kg / mm^2");
            Console.WriteLine("усилие на нижней матр = " + pn + " kg");

            Console.WriteLine("Радиал. матрица без прижима 2матр");
            Console.WriteLine("Радиус матрицы = " + rm);
        }
        //На вход: диаметры пуансонов, диаметры нижних матриц, коэфф. утон верхней матр, коэфф. утон. нижней матр, суммарный коэфф. утон на операции
        //На выход: Коэф утонений   
        public static void koeffUton2матр()
        {
            int n = int.Parse(Console.ReadLine());

            // Create arrays to store the diameters and coefficients
            int[] dp = new int[n - 1];
            int[] dm = new int[n - 1];
            float[] msw = new float[n - 1];
            float[] msn = new float[n - 1];
            float[] ms = new float[n];

            // Read the data from the console
            for (int i = 0; i < n - 1; i++)
            {
                dp[i] = int.Parse(Console.ReadLine());
                dm[i] = int.Parse(Console.ReadLine());
            }

            // Calculate the total coefficient of attenuation for each operation
            for (int j = 0; j < n; j++)
            {
                ms[j] = msw[j] * msn[j];
            }
        }
        //На вход: диамПредыдущейВыт, диамМатр, радиусМатр
        //На выход: уголРадиалМатр
        public double AngleOfRadialMatrix(double previousDiameter, double diameter, double radius)
        {
            double d1 = previousDiameter;
            double d2 = diameter;
            double r = radius;
            double a = r - (d1 - d2) / 2;
            double b = Math.Sqrt(r * r - a * a);
            double c = (d1 - d2) / 2;
            double alfa = Math.Atan(c / b);
            alfa *= 180 / 3.1416;
            return alfa;
        }
        //На вход: коэфВытяжки, уголКонусаМатр, коэфТрения, исходнаяТолщЛенты
        //На выход: толщСтенкиКолпУДна
        public double ThicknessOfConeWall(double drawRatio, double coneAngle, double frictionCoefficient, double initialThickness)
        {
            double md = drawRatio;
            double alfa = coneAngle;
            double Pi = 3.1416;
            alfa *= Pi / 180;
            double mu = frictionCoefficient;
            double s0 = initialThickness;
            double a = 1 + mu / Math.Tan(alfa);
            double b = 1 / md - 1;
            double c = 0.5 - 0.75 * a * b / (2 - 0.5 * a * b);
            double s = s0 * Math.Pow(1 / md, c);
            return s;
        }
        //На вход: Предельный коэф. утонения для стали и алюмю=0,4 для латуни=0,35, угол конусности матрицы, dn, deltan, dw. deltaw
        //На выход: средний диаметр свёртки, толщина стенки свёртки, исходная толщина равная толщине дна
        public void ВырВытяжкаСдвигом()
        {
            double dn = double.Parse(Console.ReadLine());
            double deltan = double.Parse(Console.ReadLine());
            double dw = double.Parse(Console.ReadLine());
            double deltaw = double.Parse(Console.ReadLine());
            double d1 = (dn + deltan / 2 + dw + deltaw) / 2;
            double s1 = ((dn + deltan / 2) - (dw + deltaw));
            double s0 = double.Parse(Console.ReadLine());
            double sd1 = s1 / d1;
            Console.WriteLine("Предельный коэф. утонения для стали и алюмю=0,4 для латуни=0,35");
            Console.WriteLine("Введи угол конусности матрицы");
            double alf = double.Parse(Console.ReadLine());
            Console.WriteLine(sd1);
            double kt1 = s0 / s1;
            double al = 3.1416 * alf / 180;
            double md1 = 1 / (1 + sd1 * (2 - 1 / kt1));
            double m = 1 + (1 - md1 * md1) * kt1 * kt1 * Math.Tan(al) / (2 * sd1 * md1 * md1);
            double ms1 = 1 / Math.Sqrt(m);
            Console.WriteLine(ms1);
            Console.WriteLine("штамп. сдвигом");
        }
    }
}


