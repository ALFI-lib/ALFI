#include <qcustomplot.h>

#include <ALFI.h>

template <typename T>
QVector<T> to_qvector(const std::vector<T>& v) {
	return QVector<T>(v.begin(), v.end());
}

double f(double x) {
	// return -3 * std::sin(10 * x) + std::abs(x) + 0.5 * x - x * x;
	return -3 * std::sin(10 * x) + 10 * sin(std::abs(x) + 0.5 * x);
}

template <typename Container>
Container f(const Container& X) {
	Container Y;
	Y.reserve(X.size());
	for (const auto& x : X) {
		Y.push_back(f(x));
	}
	return Y;
}

class PlotWindow final : public QWidget {
public:
	PlotWindow() {
		_plot = new QCustomPlot();
		_plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

		_control_panel = new QWidget();
		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* control_layout = new QVBoxLayout(_control_panel);

		_n_spin_box = create_spin_box<int, QSpinBox>(11, 1);

		_nn_spin_box = create_spin_box<int, QSpinBox>(10000, 1);

		_a_spin_box = create_spin_box<double, QDoubleSpinBox>(-10, -0.1);
		_b_spin_box = create_spin_box<double, QDoubleSpinBox>(10, 0.1);

		_distribution_combo = new QComboBox();
		_distribution_combo->addItem("Uniform");
		_distribution_combo->addItem("Chebyshev");
		_distribution_combo->addItem("Chebyshev Stretched");
		_distribution_combo->addItem("Chebyshev Ellipse");
		_distribution_combo->addItem("Chebyshev Ellipse Stretched");
		_distribution_combo->addItem("Circle Projection");
		_distribution_combo->addItem("Ellipse Projection");
		_distribution_combo->addItem("Sigmoid");
		_distribution_combo->addItem("Stretched Sigmoid");
		_distribution_combo->addItem("Error Function");
		_distribution_combo->addItem("Stretched Error Function");

		_function_checkbox = new QCheckBox("Function");
		_points_checkbox = new QCheckBox("Points");
		_lagrange_checkbox = new QCheckBox("Lagrange Polynomial");
		_newton_checkbox = new QCheckBox("Newton Polynomial");
		_barycentric_checkbox = new QCheckBox("Barycentric Formula");

		_function_checkbox->setChecked(true);
		_points_checkbox->setChecked(true);
		_barycentric_checkbox->setChecked(true);

		_distribution_barycentric_combo = new QComboBox();
		_distribution_barycentric_combo->addItem("Auto");
		_distribution_barycentric_combo->addItem("Uniform");
		_distribution_barycentric_combo->addItem("Chebyshev");
		_distribution_barycentric_combo->addItem("Chebyshev Stretched");
		_distribution_barycentric_combo->addItem("Chebyshev Ellipse");
		_distribution_barycentric_combo->addItem("Chebyshev Ellipse Stretched");
		_distribution_barycentric_combo->addItem("Circle Projection");
		_distribution_barycentric_combo->addItem("Ellipse Projection");
		_distribution_barycentric_combo->addItem("Sigmoid");
		_distribution_barycentric_combo->addItem("Stretched Sigmoid");
		_distribution_barycentric_combo->addItem("Error Function");
		_distribution_barycentric_combo->addItem("Stretched Error Function");

		_view_reset_button = new QPushButton("Reset View");

		connect(_n_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);

		connect(_nn_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);

		connect(_a_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_b_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);

		connect(_function_checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
		connect(_points_checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
		connect(_lagrange_checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
		connect(_newton_checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
		connect(_barycentric_checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);

		connect(_distribution_combo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &PlotWindow::update_plot);
		connect(_distribution_barycentric_combo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &PlotWindow::update_plot);

		connect(_view_reset_button, &QPushButton::clicked, this, [this]() {
			_default_axis_ranges = true;
			update_plot();
		});

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* n_layout = new QVBoxLayout();
		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* n_label_spin_box_layout = new QHBoxLayout();
		n_label_spin_box_layout->addWidget(new QLabel("N:"));
		n_label_spin_box_layout->addWidget(_n_spin_box);
		n_layout->addLayout(n_label_spin_box_layout);
		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* nn_label_spin_box_layout = new QHBoxLayout();
		nn_label_spin_box_layout->addWidget(new QLabel("nn:"));
		nn_label_spin_box_layout->addWidget(_nn_spin_box);
		n_layout->addLayout(nn_label_spin_box_layout);

		// ReSharper disable once CppDFAMemoryLeak
		QHBoxLayout* a_b_layout = new QHBoxLayout();
		a_b_layout->addWidget(new QLabel("a:"));
		a_b_layout->addWidget(_a_spin_box);
		a_b_layout->addWidget(new QLabel("b:"));
		a_b_layout->addWidget(_b_spin_box);

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* interpolation_layout = new QVBoxLayout();
		interpolation_layout->addWidget(_function_checkbox);
		interpolation_layout->addWidget(_points_checkbox);
		interpolation_layout->addWidget(_lagrange_checkbox);
		interpolation_layout->addWidget(_newton_checkbox);
		interpolation_layout->addWidget(_barycentric_checkbox);
		interpolation_layout->addWidget(_distribution_barycentric_combo);

		control_layout->addLayout(n_layout);
		control_layout->addLayout(a_b_layout);
		control_layout->addLayout(interpolation_layout);
		control_layout->addWidget(_distribution_combo);
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
	static SpinBoxType* create_spin_box(Number initial, Number step) {
		SpinBoxType* spin_box = new SpinBoxType();
		spin_box->setRange(std::min(std::numeric_limits<Number>::min(), -std::numeric_limits<Number>::infinity()),
			std::max(std::numeric_limits<Number>::max(), std::numeric_limits<Number>::infinity()));
		spin_box->setValue(initial);
		spin_box->setSingleStep(step);
		return spin_box;
	}

	void update_plot() {
		const size_t N = _n_spin_box->value();
		const size_t nn = _nn_spin_box->value();
		const double a = _a_spin_box->value();
		const double b = _b_spin_box->value();
		const int dist_type = _distribution_combo->currentIndex() + 1; // + 1 because GENERAL is first

		std::vector<double> X;
		switch (static_cast<alfi::dist::Type>(dist_type)) {
		case alfi::dist::Type::UNIFORM: X = alfi::dist::uniform(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV: X = alfi::dist::chebyshev(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV_STRETCHED: X = alfi::dist::chebyshev_stretched(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV_ELLIPSE: X = alfi::dist::chebyshev_ellipse(N, a, b, 2.0); break;
			case alfi::dist::Type::CHEBYSHEV_ELLIPSE_STRETCHED: X = alfi::dist::chebyshev_ellipse_stretched(N, a, b, 2.0); break;
			case alfi::dist::Type::CIRCLE_PROJECTION: X = alfi::dist::circle_proj(N, a, b); break;
			case alfi::dist::Type::ELLIPSE_PROJECTION: X = alfi::dist::ellipse_proj(N, a, b, 2.0); break;
			case alfi::dist::Type::SIGMOID: X = alfi::dist::sigmoid(N, a, b, 16.0); break;
			case alfi::dist::Type::SIGMOID_STRETCHED: X = alfi::dist::sigmoid_stretched(N, a, b, 16.0); break;
			case alfi::dist::Type::ERF: X = alfi::dist::erf(N, a, b, 8.0); break;
			case alfi::dist::Type::ERF_STRETCHED: X = alfi::dist::erf_stretched(N, a, b, 8.0); break;
			default: return;
		}

		const std::vector<double> Y = f(X);

		const std::vector<double> xx = alfi::dist::uniform(nn, a, b);

		double y_min = *std::ranges::min_element(Y);
		double y_max = *std::ranges::max_element(Y);

		_plot->clearGraphs();
		int graph_index = -1;

		if (_function_checkbox->isChecked()) {
			_plot->addGraph();
			++graph_index;
			const auto yy = f(xx);
			y_min = std::min(y_min, *std::ranges::min_element(yy));
			y_max = std::max(y_max, *std::ranges::max_element(yy));
			_plot->graph(graph_index)->setData(to_qvector(xx), to_qvector(yy));
			_plot->graph(graph_index)->setPen(QPen(Qt::blue));
		}

		if (_points_checkbox->isChecked()) {
			_plot->addGraph();
			++graph_index;
			_plot->graph(graph_index)->setData(to_qvector(X), to_qvector(Y));
			_plot->graph(graph_index)->setLineStyle(QCPGraph::lsNone);
			_plot->graph(graph_index)->setScatterStyle(QCPScatterStyle::ssCircle);
			_plot->graph(graph_index)->setPen(QPen(Qt::red));
		}

		if (_lagrange_checkbox->isChecked()) {
			_plot->addGraph();
			graph_index++;
			const auto coeffs = alfi::poly::lagrange(X, Y);
			const auto yy = alfi::poly::val(coeffs, xx);
			_plot->graph(graph_index)->setData(to_qvector(xx), to_qvector(yy));
			_plot->graph(graph_index)->setPen(QPen(Qt::green));
		}

		if (_newton_checkbox->isChecked()) {
			_plot->addGraph();
			graph_index++;
			const auto coeffs = alfi::poly::newton(X, Y);
			const auto yy = alfi::poly::val(coeffs, xx);
			_plot->graph(graph_index)->setData(to_qvector(xx), to_qvector(yy));
			_plot->graph(graph_index)->setPen(QPen(Qt::red));
		}

		if (_barycentric_checkbox->isChecked()) {
			_plot->addGraph();
			graph_index++;
			int barycentric_dist_type = _distribution_barycentric_combo->currentIndex();
			if (barycentric_dist_type == 0) {
				barycentric_dist_type = dist_type;
			}
			const auto yy = alfi::misc::barycentric(X, Y, xx, static_cast<alfi::dist::Type>(barycentric_dist_type));
			_plot->graph(graph_index)->setData(to_qvector(xx), to_qvector(yy));
			_plot->graph(graph_index)->setPen(QPen(Qt::magenta));
		}

		if (_default_axis_ranges) {
			const double RESERVE = 0.1;
			const double FACTOR = 1 + RESERVE;

			const double x_center = (a + b) / 2;
			const double x_radius = (b - a) / 2;
			const double x_min = x_center - FACTOR * x_radius;
			const double x_max = x_center + FACTOR * x_radius;
			_plot->xAxis->setRange(x_min, x_max);

			const double y_center = (y_min + y_max) / 2;
			const double y_radius = (y_max - y_min) / 2;
			const double y_min_with_padding = y_center - FACTOR * y_radius;
			const double y_max_with_padding = y_center + FACTOR * y_radius;
			_plot->yAxis->setRange(y_min_with_padding, y_max_with_padding);
			_default_axis_ranges = false;
		}

		_plot->replot();
	}

private:
	bool _default_axis_ranges = true;
	QCustomPlot* _plot;
	QWidget* _control_panel;
	QSpinBox* _n_spin_box;
	QSpinBox* _nn_spin_box;
	QDoubleSpinBox* _a_spin_box;
	QDoubleSpinBox* _b_spin_box;
	QCheckBox* _function_checkbox;
	QCheckBox* _points_checkbox;
	QCheckBox* _lagrange_checkbox;
	QCheckBox* _newton_checkbox;
	QCheckBox* _barycentric_checkbox;
	QComboBox* _distribution_barycentric_combo;
	QComboBox* _distribution_combo;
	QPushButton* _view_reset_button;
};

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);

	PlotWindow window;
	window.setWindowTitle("Interpolation");
	window.resize(1200, 600);
	window.show();

	return QApplication::exec();
}