using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace Mothra.UI
{
    /// <summary>
    /// FigureUI.xaml の相互作用ロジック
    /// </summary>
    public partial class ControlBox : Window
    {
        public System.Collections.Generic.List<Mothra.UI.slider> listSlider;
        public newRadioButton[] radioButtonList;
        public void setNumF(int num)
        {
            radioButtonList=new newRadioButton[num];
            for (int i = 0; i < num; i++)
            {
                radioButtonList[i] = new newRadioButton();
                radioButtonList[i].number = i;
                this.wrapPanel1.Children.Add(radioButtonList[i]);
                radioButtonList[i].radiobutton1.Content = i.ToString();
            }
            foreach (var rB in radioButtonList)
            {
                rB.radiobutton1.IsEnabled = false;
                rB.radiobutton1.GroupName = "funcs";
            }
        }
        public void setFunctionToCompute(Action func)
        {
            this.Compute.function = func;
        }
        public void setFunctionToCompute2(Action func)
        {
            this.Compute2.function = func;
        }
        public ControlBox()
        {
            InitializeComponent();
            listSlider = new System.Collections.Generic.List<Mothra.UI.slider>();
            this.Hide();
        }
        public void clearSliders()
        {
            this.stackPanel1.Children.Clear();
            listSlider.Clear();
        }
        public Mothra.UI.newButton addButton(string name,Action func)
        {
            Mothra.UI.newButton b=new Mothra.UI.newButton(func);
            this.stackPanel1.Children.Add(b);
            b.Text = name;
            return b;
        }
        public Mothra.UI.slider addSlider(int min, int step, int max,int val)
        {
            Mothra.UI.slider s = addSlider(min, step, max, val, "");
            return s;
        }
        public Mothra.UI.slider addSlider(int min, int step, int max, int val, string text)
        {
            Mothra.UI.slider s = new Mothra.UI.slider(min, step, max, val, text);
            this.stackPanel1.Children.Add(s);
            s.getSlider.ValueChanged +=new RoutedPropertyChangedEventHandler<double>(this.slider_ValueChanged);
            listSlider.Add(s);
            return s;
        }

        private void slider_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            ((Mothra.UI.slider)(from p in listSlider where (System.Object.ReferenceEquals(p.getSlider,sender))select p).First()).update();
        }

        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            e.Cancel = true;
            this.Hide();
        }

        public void EnableRadio(int lastComputed,Action<int> func)
        {
            radioButtonList[lastComputed].radiobutton1.IsEnabled = true;
            radioButtonList[lastComputed].function = func;
            radioButtonList[lastComputed].number = lastComputed;
        }
    }
}
