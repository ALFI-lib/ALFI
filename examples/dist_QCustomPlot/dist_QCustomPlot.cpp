#include <qcustomplot.h>

#include <ALFI.h>

class PlotWindow final : public QWidget {
private:
	template <typename Number>
	struct VariableParams {
		Number initial;
		Number min;
		Number max;
		Number step;
	};

public:
	PlotWindow() {
		static const VariableParams<int> n_params = {7, 1, 40, 1};

		static const VariableParams<double> a_params = {0.0, 0.0, 100.0, 0.01};
		static const VariableParams<double> b_params = {1.0, 0.0, 100.0, 0.01};
		static const VariableParams<double> ratio_params = {2.0, 0.0, 4.0, 0.01};
		static const VariableParams<double> steepness_params = {8.0, 0.0, 10.0, 0.01};

		_plot = new QCustomPlot();
		_plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
		_plot->legend->setVisible(true);
		_plot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));

		_control_panel = new QWidget();
		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* control_layout = new QVBoxLayout(_control_panel);

		_n_spin_box = create_spin_box<int, QSpinBox>(n_params);
		_n_slider = create_slider<int, QSlider>(n_params);

		_a_spin_box = create_spin_box<double, QDoubleSpinBox>(a_params);
		_b_spin_box = create_spin_box<double, QDoubleSpinBox>(b_params);

		_ratio_spin_box = create_spin_box<double, QDoubleSpinBox>(ratio_params);
		_ratio_slider = create_slider<double, QSlider>(ratio_params);

		_steepness_spin_box = create_spin_box<double, QDoubleSpinBox>(steepness_params);
		_steepness_slider = create_slider<double, QSlider>(steepness_params);

		_view_reset_button = new QPushButton("Reset View");

		if (!connect(_n_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int value) {
			_n_slider->setValue(value);
			update_plot();
		}))
			qDebug() << "Failed to connect _n_spin_box valueChanged signal";
		if (!connect(_a_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot))
			qDebug() << "Failed to connect _a_spin_box valueChanged signal";
		if (!connect(_b_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot))
			qDebug() << "Failed to connect _b_spin_box valueChanged signal";
		if (!connect(_ratio_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double value) {
			_ratio_slider->setValue(static_cast<int>(value / ratio_params.step));
			update_plot();
		}))
			qDebug() << "Failed to connect _ratio_spin_box valueChanged signal";
		if (!connect(_steepness_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double value) {
			_steepness_slider->setValue(static_cast<int>(value / steepness_params.step));
			update_plot();
		}))
			qDebug() << "Failed to connect _steepness_spin_box valueChanged signal";

		if (!connect(_n_slider, &QSlider::valueChanged, this, [this](int value) {
			_n_spin_box->setValue(value);
			update_plot();
		}))
			qDebug() << "Failed to connect _n_slider valueChanged signal";
		if (!connect(_ratio_slider, &QSlider::valueChanged, this, [this](int value) {
			_ratio_spin_box->setValue(value * ratio_params.step);
			update_plot();
		}))
			qDebug() << "Failed to connect _ratio_slider valueChanged signal";
		if (!connect(_steepness_slider, &QSlider::valueChanged, this, [this](int value) {
			_steepness_spin_box->setValue(value * steepness_params.step);
			update_plot();
		}))
			qDebug() << "Failed to connect _steepness_slider valueChanged signal";
		if (!connect(_view_reset_button, &QPushButton::clicked, this, [this]() {
			_default_axis_ranges = true;
			update_plot();
		}))
			qDebug() << "Failed to connect _view_reset_button clicked signal";

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* n_layout = new QVBoxLayout();
		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* n_label_spin_box_layout = new QHBoxLayout();
		n_label_spin_box_layout->addWidget(new QLabel("n:"));
		n_label_spin_box_layout->addWidget(_n_spin_box);
		n_layout->addLayout(n_label_spin_box_layout);
		n_layout->addWidget(_n_slider);

		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* a_b_layout = new QHBoxLayout();
		a_b_layout->addWidget(new QLabel("a:"));
		a_b_layout->addWidget(_a_spin_box);
		a_b_layout->addWidget(new QLabel("b:"));
		a_b_layout->addWidget(_b_spin_box);

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* B_layout = new QVBoxLayout();
		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* B_label_spin_box_layout = new QHBoxLayout();
		B_label_spin_box_layout->addWidget(new QLabel("Ellipse ratio:"));
		B_label_spin_box_layout->addWidget(_ratio_spin_box);
		B_layout->addLayout(B_label_spin_box_layout);
		B_layout->addWidget(_ratio_slider);

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* steepness_layout = new QVBoxLayout();
		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* steepness_label_spin_box_layout = new QHBoxLayout();
		steepness_label_spin_box_layout->addWidget(new QLabel("Steepness:"));
		steepness_label_spin_box_layout->addWidget(_steepness_spin_box);
		steepness_layout->addLayout(steepness_label_spin_box_layout);
		steepness_layout->addWidget(_steepness_slider);

		control_layout->addLayout(n_layout);
		control_layout->addLayout(a_b_layout);
		control_layout->addLayout(B_layout);
		control_layout->addLayout(steepness_layout);
		control_layout->addWidget(_view_reset_button);

		_control_panel->setMinimumWidth(250);

		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* main_layout = new QHBoxLayout();
		main_layout->addWidget(_plot, 1);
		main_layout->addWidget(_control_panel);

		setLayout(main_layout);

		update_plot();
	}

private:
	template <typename Number, typename SpinBoxType>
	static SpinBoxType* create_spin_box(const VariableParams<Number>& params) {
		SpinBoxType* spin_box = new SpinBoxType();
		spin_box->setRange(params.min, params.max);
		spin_box->setValue(params.initial);
		spin_box->setSingleStep(params.step);
		return spin_box;
	}

	template <typename Number, typename SliderType>
	static SliderType* create_slider(const VariableParams<Number>& params) {
		SliderType* slider = new SliderType(Qt::Horizontal);
		slider->setRange(static_cast<int>(params.min / params.step), static_cast<int>(params.max / params.step));
		slider->setValue(static_cast<int>(params.initial / params.step));
		return slider;
	}

	void update_plot() {
		static const QVector<QColor> colors = {
			QColor(  0,   0, 255), // blue        #0000ff =   0   0 255
			QColor(255,   0,   0), // red         #ff0000 = 255   0   0
			QColor(  0, 255,   0), // green       #00ff00 =   0 255   0
			QColor(255, 165,   0), // orange      #ffa500 = 255 165   0
			QColor(192, 128, 255), // purple      #c080ff = 192 128 255
			QColor(165,  42,  42), // brown       #a52a2a = 165  42  42
			QColor(255,   0, 255), // magenta     #ff00ff = 255   0 255
			QColor(  0,   0, 139), // dark-blue   #00008b =   0   0 139
			QColor(255, 127,  80), // coral       #ff7f50 = 255 127  80
			QColor(148,   0, 211), // dark-violet #9400d3 = 148   0 211
			QColor(250, 128, 114), // salmon      #fa8072 = 250 128 114
			QColor(139,   0,   0), // dark-red    #8b0000 = 139   0   0
			QColor(  0, 128, 255), // web-blue    #0080ff =   0 128 255
			QColor(  0, 192,   0), // web-green   #00c000 =   0 192   0
			QColor(192,  64,   0), // dark-orange #c04000 = 192  64   0
		};

		const size_t n = _n_spin_box->value();
		const double a = _a_spin_box->value();
		const double b = _b_spin_box->value();
		const double B = _ratio_spin_box->value();
		const double steepness = _steepness_spin_box->value();

		const QVector<std::pair<QString, std::vector<double>>> distributions = {
			{"Uniform", alfi::dist::uniform(n, a, b)},
			{"Quadratic", alfi::dist::quadratic(n, a, b)},
			{"Cubic", alfi::dist::cubic(n, a, b)},
			{"Chebyshev", alfi::dist::chebyshev(n, a, b)},
			{"Stretched Chebyshev", alfi::dist::chebyshev_stretched(n, a, b)},
			{"Circle Projection", alfi::dist::circle_proj(n, a, b)},
			{"Chebyshev Ellipse", alfi::dist::chebyshev_ellipse(n, a, b, B)},
			{"Stretched Chebyshev Ellipse", alfi::dist::chebyshev_ellipse_stretched(n, a, b, B)},
			{"Ellipse Projection", alfi::dist::ellipse_proj(n, a, b, B)},
			{"Logistic", alfi::dist::logistic(n, a, b, steepness)},
			{"Stretched Logistic", alfi::dist::logistic_stretched(n, a, b, steepness)},
			{"Error function", alfi::dist::erf(n, a, b, steepness)},
			{"Stretched Error function", alfi::dist::erf_stretched(n, a, b, steepness)},
		};

		const QCPRange current_x_range = _plot->xAxis->range();
		const QCPRange current_y_range = _plot->yAxis->range();

		_plot->clearGraphs();

		double current_layer = -1;

		for (int i = 0; i < distributions.size(); ++i) {
			const auto& [name, x_data] = distributions[i];
			QVector<double> x, y;
			for (const double& value : x_data) {
				if (!std::isnan(value)) {
					x.append(value);
					y.append(current_layer);
				}
			}

			_plot->addGraph();
			_plot->graph(i)->setName(name);
			_plot->graph(i)->setData(x, y);
			_plot->graph(i)->setLineStyle(QCPGraph::lsNone);
			const QColor& color = colors[i % colors.size()];
			_plot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, color, color, 12));

			current_layer -= 1;
		}

		if (_default_axis_ranges) {
			_plot->xAxis->setRange(a, b);
			_plot->yAxis->setRange(current_layer, 0);
			_default_axis_ranges = false;
		} else {
			_plot->xAxis->setRange(current_x_range);
			_plot->yAxis->setRange(current_y_range);
		}

		_plot->replot();
	}

private:
	bool _default_axis_ranges = true;
	QCustomPlot* _plot;
	QWidget* _control_panel;
	QSpinBox* _n_spin_box;
	QSlider* _n_slider;
	QDoubleSpinBox* _a_spin_box;
	QDoubleSpinBox* _b_spin_box;
	QDoubleSpinBox* _ratio_spin_box;
	QSlider* _ratio_slider;
	QDoubleSpinBox* _steepness_spin_box;
	QSlider* _steepness_slider;
	QPushButton* _view_reset_button;
};

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);

	PlotWindow window;
	window.setWindowTitle("ALFI Plot with Controls");
	window.resize(1200, 600);
	window.show();

	return QApplication::exec();
}