#include <ranges>

#include <qcustomplot.h>

#include <ALFI.h>

template <typename T>
QVector<T> to_qvector(const std::vector<T>& v) {
	return QVector<T>(v.begin(), v.end());
}

double f(double x) {
	// return -3 * std::sin(10 * x) + std::abs(x) + 0.5 * x - x * x; // example 1
	return -3 * std::sin(10 * x) + 10 * sin(std::abs(x) + 0.5 * x); // example 2
	// return std::sin(3 * M_PI * x) + std::cos(M_PI * x / 5); // periodic
	// return 1/(1+std::exp(-10000*(x-1))) - 1/(1+std::exp(-10000*(x+1))); // high local second derivative
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
	static const inline std::vector<std::pair<QString,std::function<std::vector<double>(std::vector<double>,std::vector<double>,std::vector<double>)>>>
	poly_types = {
		{"Lagrange", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::val(alfi::poly::lagrange(X, Y), xx); }},
		{"Lagrange values", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::lagrange_vals(X, Y, xx); }},
		{"Improved Lagrange", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::val(alfi::poly::imp_lagrange(X, Y), xx); }},
		{"Improved Lagrange values", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::imp_lagrange_vals(X, Y, xx); }},
		{"Newton", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::val(alfi::poly::newton(X, Y), xx); }},
		{"Newton values", [](const auto& X, const auto& Y, const auto& xx) { return alfi::poly::newton_vals(X, Y, xx); }},
	};
	static const inline size_t poly_default_type_index = 0;

	static const inline std::vector<std::pair<QString,alfi::spline::StepSpline<>::Type>> step_spline_types = {
		{"Left", alfi::spline::StepSpline<>::Types::Left{}},
		{"Middle", alfi::spline::StepSpline<>::Types::Middle{}},
		{"Right", alfi::spline::StepSpline<>::Types::Right{}},
	};
	static const inline size_t step_spline_default_type_index = 0;

	static const inline std::vector<std::pair<QString,alfi::spline::QuadraticSpline<>::Type>> quadratic_spline_types = {
		{"Not-a-knot start", alfi::spline::QuadraticSpline<>::Types::NotAKnotStart{}},
		{"Not-a-knot end", alfi::spline::QuadraticSpline<>::Types::NotAKnotEnd{}},
		{"Semi-not-a-knot", alfi::spline::QuadraticSpline<>::Types::SemiNotAKnot{}},
		{"Natural start", alfi::spline::QuadraticSpline<>::Types::NaturalStart{}},
		{"Natural end", alfi::spline::QuadraticSpline<>::Types::NaturalEnd{}},
		{"Semi-natural", alfi::spline::QuadraticSpline<>::Types::SemiNatural{}},
		{"Semi-semi", alfi::spline::QuadraticSpline<>::Types::SemiSemi{}},
		{"Clamped-start(10)", alfi::spline::QuadraticSpline<>::Types::ClampedStart{10}},
		{"Clamped-start(-10)", alfi::spline::QuadraticSpline<>::Types::ClampedStart{-10}},
		{"Clamped-end(10)", alfi::spline::QuadraticSpline<>::Types::ClampedEnd{10}},
		{"Clamped-end(-10)", alfi::spline::QuadraticSpline<>::Types::ClampedEnd{-10}},
		{"Semi-clamped(10, 10)", alfi::spline::QuadraticSpline<>::Types::SemiClamped{10, 10}},
		{"Fixed-second-start(10)", alfi::spline::QuadraticSpline<>::Types::FixedSecondStart{10}},
		{"Fixed-second-end(10)", alfi::spline::QuadraticSpline<>::Types::FixedSecondEnd{10}},
		{"Semi-fixed-second(10, 10)", alfi::spline::QuadraticSpline<>::Types::SemiFixedSecond{10, 10}},
		{"Clamped(10, 10)", alfi::spline::QuadraticSpline<>::Types::Clamped{10, 10}},
		{"Fixed-second(10, 10)", alfi::spline::QuadraticSpline<>::Types::FixedSecond{10, 10}},
		{"Not-a-knot(10)", alfi::spline::QuadraticSpline<>::Types::NotAKnot{10}},
	};
	static const inline size_t quadratic_spline_default_type_index = 2;

	static const inline std::vector<std::pair<QString,alfi::spline::CubicSpline<>::Type>> cubic_spline_types = {
		{"Natural", alfi::spline::CubicSpline<>::Types::Natural{}},
		{"Not-a-knot", alfi::spline::CubicSpline<>::Types::NotAKnot{}},
		{"Periodic", alfi::spline::CubicSpline<>::Types::Periodic{}},
		{"ParabolicEnds", alfi::spline::CubicSpline<>::Types::ParabolicEnds{}},
		{"Clamped(0, 10)", alfi::spline::CubicSpline<>::Types::Clamped{0, 10}},
		{"FixedSecond(40, 20)", alfi::spline::CubicSpline<>::Types::FixedSecond{40, 20}},
		{"FixedThird(10, 50)", alfi::spline::CubicSpline<>::Types::FixedThird{10, 50}},
		{"Not-a-knot start", alfi::spline::CubicSpline<>::Types::NotAKnotStart{}},
		{"Not-a-knot end", alfi::spline::CubicSpline<>::Types::NotAKnotEnd{}},
		{"Semi-not-a-knot", alfi::spline::CubicSpline<>::Types::SemiNotAKnot{}},
	};
	static const inline size_t cubic_spline_default_type_index = 0;

	PlotWindow() {
		static const QStringList distribution_types {
			"Uniform", "Quadratic", "Cubic", "Chebyshev", "Stretched Chebyshev",
			"Circle Projection", "Circle Projection Without Last", "Circle Projection Without First",
			"Chebyshev Ellipse", "Stretched Chebyshev Ellipse", "Ellipse Projection",
			"Ellipse Projection Without Last", "Ellipse Projection Without First",
			"Logistic", "Stretched Logistic", "Error Function", "Stretched Error Function"
		};

		_plot = new QCustomPlot();
		_plot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
		_plot->setPlottingHints(QCP::phFastPolylines | QCP::phCacheLabels);
		_plot->setAntialiasedElements(QCP::aeAll);
		_plot->legend->setVisible(true);
		_plot->legend->setBrush(QBrush(QColor(255, 255, 255, 0)));

		_control_panel = new QWidget();
		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* control_layout = new QVBoxLayout(_control_panel);

		_n_spin_box = create_spin_box<int, QSpinBox>(11, 1, 1);

		_nn_spin_box = create_spin_box<int, QSpinBox>(10000, 1, 1);

		_a_spin_box = create_spin_box<double, QDoubleSpinBox>(-10, -0.1);
		_b_spin_box = create_spin_box<double, QDoubleSpinBox>(10, 0.1);

		_distribution_combo = new QComboBox();
		_distribution_combo->addItems(distribution_types);

		_function_checkbox = new QCheckBox("Function");
		_points_checkbox = new QCheckBox("Points");
		_poly_checkbox = new QCheckBox("Polynomial");
		_barycentric_checkbox = new QCheckBox("Barycentric Formula");
		_poly_eqv_spline_checkbox = new QCheckBox("Poly. Eqv. Spline");
		_step_spline_checkbox = new QCheckBox("Step Spline");
		_linear_spline_checkbox = new QCheckBox("Linear Spline");
		_quadratic_spline_checkbox = new QCheckBox("Quadratic Spline");
		_cubic_spline_checkbox = new QCheckBox("Cubic Spline");

		_function_checkbox->setChecked(true);
		_points_checkbox->setChecked(true);

		_poly_combo = new QComboBox();
		for (const auto& name : poly_types | std::views::keys) {
			_poly_combo->addItem(name);
		}
		_poly_combo->setCurrentIndex(poly_default_type_index);

		_barycentric_combo = new QComboBox();
		_barycentric_combo->addItem("Auto");
		_barycentric_combo->addItems(distribution_types);

		_step_spline_combo = new QComboBox();
		for (const auto& name : step_spline_types | std::views::keys) {
			_step_spline_combo->addItem(name);
		}
		_step_spline_combo->setCurrentIndex(step_spline_default_type_index);

		_quadratic_spline_combo = new QComboBox();
		for (const auto& name : quadratic_spline_types | std::views::keys) {
			_quadratic_spline_combo->addItem(name);
		}
		_quadratic_spline_combo->setCurrentIndex(quadratic_spline_default_type_index);

		_cubic_spline_combo = new QComboBox();
		for (const auto& name : cubic_spline_types | std::views::keys) {
			_cubic_spline_combo->addItem(name);
		}
		_cubic_spline_combo->setCurrentIndex(cubic_spline_default_type_index);

		_view_reset_button = new QPushButton("Reset View");

		connect(_n_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_nn_spin_box, QOverload<int>::of(&QSpinBox::valueChanged), this, &PlotWindow::update_plot);

		connect(_a_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);
		connect(_b_spin_box, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &PlotWindow::update_plot);

		for (const auto checkbox : {
			_function_checkbox, _points_checkbox, _poly_checkbox, _barycentric_checkbox, _poly_eqv_spline_checkbox,
			_step_spline_checkbox, _linear_spline_checkbox, _quadratic_spline_checkbox, _cubic_spline_checkbox}) {
			connect(checkbox, &QCheckBox::toggled, this, &PlotWindow::update_plot);
		}

		for (const auto combo : {_distribution_combo, _poly_combo, _barycentric_combo, _step_spline_combo, _quadratic_spline_combo, _cubic_spline_combo}) {
			connect(combo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &PlotWindow::update_plot);
		}

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
		// ReSharper disable once CppDFAMemoryLeak
		QLabel* a_label = new QLabel("a:");
		a_label->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
		_a_spin_box->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
		// ReSharper disable once CppDFAMemoryLeak
		QLabel* b_label = new QLabel("b:");
		b_label->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
		_b_spin_box->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
		a_b_layout->addWidget(a_label);
		a_b_layout->addWidget(_a_spin_box);
		a_b_layout->addWidget(b_label);
		a_b_layout->addWidget(_b_spin_box);

		// ReSharper disable once CppDFAMemoryLeak
		QVBoxLayout* interpolation_layout = new QVBoxLayout();
		interpolation_layout->addWidget(_function_checkbox);
		interpolation_layout->addWidget(_points_checkbox);
		interpolation_layout->addWidget(_poly_checkbox);
		interpolation_layout->addWidget(_poly_combo);
		interpolation_layout->addWidget(_barycentric_checkbox);
		interpolation_layout->addWidget(_barycentric_combo);
		interpolation_layout->addWidget(_poly_eqv_spline_checkbox);
		interpolation_layout->addWidget(_step_spline_checkbox);
		interpolation_layout->addWidget(_step_spline_combo);
		interpolation_layout->addWidget(_linear_spline_checkbox);
		interpolation_layout->addWidget(_quadratic_spline_checkbox);
		interpolation_layout->addWidget(_quadratic_spline_combo);
		interpolation_layout->addWidget(_cubic_spline_checkbox);
		interpolation_layout->addWidget(_cubic_spline_combo);

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
	static SpinBoxType* create_spin_box(
			Number initial,
			Number step,
			Number min = std::min(std::numeric_limits<Number>::lowest(), -std::numeric_limits<Number>::infinity()),
			Number max = std::max(std::numeric_limits<Number>::max(), std::numeric_limits<Number>::infinity())) {
		SpinBoxType* spin_box = new SpinBoxType();
		spin_box->setRange(min, max);
		spin_box->setValue(initial);
		spin_box->setSingleStep(step);
		return spin_box;
	}

	void update_plot() {
		static const std::vector<QColor> colors = {
			QColor(  0, 102, 204), // deep bue
			QColor(204,   0,   0), // red
			QColor(  0, 153,   0), // green
			QColor(204, 102,   0), // muted orange
			QColor(102,   0, 153), // purple
			QColor( 30,  30,  30), // darker gray
			QColor(  0, 128, 128), // dark turquoise
			QColor(204,   0, 102), // pink
		};
		static const qreal line_width = 1.75;

		const size_t N = _n_spin_box->value();
		const size_t nn = _nn_spin_box->value();
		const double a = _a_spin_box->value();
		const double b = _b_spin_box->value();
		const int dist_type = _distribution_combo->currentIndex() + 1; // + 1 because GENERAL is first

		std::vector<double> X;
		switch (static_cast<alfi::dist::Type>(dist_type)) {
			case alfi::dist::Type::UNIFORM: X = alfi::dist::uniform(N, a, b); break;
			case alfi::dist::Type::QUADRATIC: X = alfi::dist::quadratic(N, a, b); break;
			case alfi::dist::Type::CUBIC: X = alfi::dist::cubic(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV: X = alfi::dist::chebyshev(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV_STRETCHED: X = alfi::dist::chebyshev_stretched(N, a, b); break;
			case alfi::dist::Type::CIRCLE_PROJ: X = alfi::dist::circle_proj(N, a, b); break;
			case alfi::dist::Type::CIRCLE_PROJ_NO_LAST: X = alfi::dist::circle_proj_no_last(N, a, b); break;
			case alfi::dist::Type::CIRCLE_PROJ_NO_FIRST: X = alfi::dist::circle_proj_no_first(N, a, b); break;
			case alfi::dist::Type::CHEBYSHEV_ELLIPSE: X = alfi::dist::chebyshev_ellipse(N, a, b, 2.0); break;
			case alfi::dist::Type::CHEBYSHEV_ELLIPSE_STRETCHED: X = alfi::dist::chebyshev_ellipse_stretched(N, a, b, 2.0); break;
			case alfi::dist::Type::ELLIPSE_PROJ: X = alfi::dist::ellipse_proj(N, a, b, 2.0); break;
			case alfi::dist::Type::ELLIPSE_PROJ_NO_LAST: X = alfi::dist::ellipse_proj_no_last(N, a, b, 2.0); break;
			case alfi::dist::Type::ELLIPSE_PROJ_NO_FIRST: X = alfi::dist::ellipse_proj_no_first(N, a, b, 2.0); break;
			case alfi::dist::Type::LOGISTIC: X = alfi::dist::logistic(N, a, b, 16.0); break;
			case alfi::dist::Type::LOGISTIC_STRETCHED: X = alfi::dist::logistic_stretched(N, a, b, 16.0); break;
			case alfi::dist::Type::ERF: X = alfi::dist::erf(N, a, b, 8.0); break;
			case alfi::dist::Type::ERF_STRETCHED: X = alfi::dist::erf_stretched(N, a, b, 8.0); break;
			default: return;
		}

		const auto Y = f(X);

		const std::vector<double> xx = alfi::dist::uniform(nn, a, b);

		double y_min = *std::ranges::min_element(Y);
		double y_max = *std::ranges::max_element(Y);

		_plot->clearGraphs();
		int graph_index = -1;

		const auto add_graph = [&](
				const auto& name,
				const std::vector<double>& x_data,
				const std::vector<double>& y_data,
				QCPGraph::LineStyle line_style = QCPGraph::lsLine,
				const QCPScatterStyle& scatter_style = QCPScatterStyle::ssNone,
				bool over = false) {
			_plot->addGraph();
			++graph_index;
			_plot->graph(graph_index)->setName(name);
			_plot->graph(graph_index)->setData(to_qvector(x_data), to_qvector(y_data));
			_plot->graph(graph_index)->setLineStyle(line_style);
			_plot->graph(graph_index)->setScatterStyle(scatter_style);
			_plot->graph(graph_index)->setPen(QPen(colors[graph_index % colors.size()], line_width));
			if (over) _plot->graph(graph_index)->setLayer("axes");
		};

		if (_function_checkbox->isChecked()) {
			const auto yy = f(xx);
			y_min = std::min(y_min, *std::ranges::min_element(yy));
			y_max = std::max(y_max, *std::ranges::max_element(yy));
			add_graph("Function", xx, yy);
		}
		if (_points_checkbox->isChecked()) {
			add_graph("Points", X, Y, QCPGraph::lsNone, QCPScatterStyle::ssCircle, true);
		}
		if (_poly_checkbox->isChecked()) {
			add_graph("Polynomial", xx, poly_types[_poly_combo->currentIndex()].second(X, Y, xx));
		}
		if (_barycentric_checkbox->isChecked()) {
			int barycentric_dist_type = _barycentric_combo->currentIndex();
			if (barycentric_dist_type == 0) {
				barycentric_dist_type = dist_type;
			}
			add_graph("Barycentric", xx, alfi::misc::barycentric(X, Y, xx, static_cast<alfi::dist::Type>(barycentric_dist_type)));
		}
		if (_poly_eqv_spline_checkbox->isChecked()) {
			add_graph("Poly. Eqv. Spline", xx, alfi::spline::PolyEqvSpline<>(X, Y)(xx));
		}
		if (_step_spline_checkbox->isChecked()) {
			add_graph("Step Spline", xx, alfi::spline::StepSpline(X, Y, step_spline_types[_step_spline_combo->currentIndex()].second)(xx));
		}
		if (_linear_spline_checkbox->isChecked()) {
			add_graph("Linear Spline", xx, alfi::spline::LinearSpline(X, Y)(xx));
		}
		if (_quadratic_spline_checkbox->isChecked()) {
			add_graph("Quadratic Spline", xx, alfi::spline::QuadraticSpline<>(X, Y, quadratic_spline_types[_quadratic_spline_combo->currentIndex()].second)(xx));
		}
		if (_cubic_spline_checkbox->isChecked()) {
			add_graph("Cubic Spline", xx, alfi::spline::CubicSpline<>(X, Y, cubic_spline_types[_cubic_spline_combo->currentIndex()].second)(xx));
		}

		if (_default_axis_ranges) {
			_default_axis_ranges = false;
			static const double RESERVE = 0.1;
			static const double FACTOR = 1 + RESERVE;
			const double x_center = (a + b) / 2, x_radius = (b - a) / 2;
			_plot->xAxis->setRange(x_center - FACTOR * x_radius, x_center + FACTOR * x_radius);
			const double y_center = (y_min + y_max) / 2, y_radius = (y_max - y_min) / 2;
			_plot->yAxis->setRange(y_center - FACTOR * y_radius, y_center + FACTOR * y_radius);
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
	QCheckBox* _poly_checkbox;
	QCheckBox* _barycentric_checkbox;
	QCheckBox* _poly_eqv_spline_checkbox;
	QCheckBox* _step_spline_checkbox;
	QCheckBox* _linear_spline_checkbox;
	QCheckBox* _quadratic_spline_checkbox;
	QCheckBox* _cubic_spline_checkbox;
	QComboBox* _distribution_combo;
	QComboBox* _poly_combo;
	QComboBox* _barycentric_combo;
	QComboBox* _step_spline_combo;
	QComboBox* _quadratic_spline_combo;
	QComboBox* _cubic_spline_combo;
	QPushButton* _view_reset_button;
};

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);

	PlotWindow window;
	window.setWindowTitle("Interpolation");
	window.resize(1280, 720);
	window.show();

	return QApplication::exec();
}