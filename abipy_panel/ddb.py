import ipywidgets as widgets
from IPython.display import display


class Summary(widgets.Output):
    def __init__(self):
        super().__init__()

    def update(self, ddb, verbose=False):
        self.clear_output()
        with self:
            print(ddb.to_string(verbose=True))


class PHBandsOptions(widgets.VBox):
    def __init__(self):
        super().__init__()

        self._header = widgets.HTML("<h1>PH-bands options</h1>")
        self._caption = widgets.HTML(
            "Compute phonon bands and DOSes from DDB and plot the results."
        )

        # "Number of divisions for smallest vector to generate Q-mesh"
        self._nqsmall = widgets.BoundedIntText(description="nqsmall:", value=10, min=1)

        # "Number of divisions for smallest vector to generate Q-path"
        self._ndivsm = widgets.BoundedIntText(
            description="ndivsm:",
            value=5,
            min=1,
        )

        self._asr = widgets.Dropdown(
            description="asr:",
            options=["0", "1", "2"],
            value="2",
        )

        self._chneut = widgets.Dropdown(
            description="chneut:",
            options=["0", "1", "2"],
            value="1",
        )

        self._dipdip = widgets.Dropdown(
            description="dipdip:",
            options=["0", "1", "2"],
            value="1",
        )

        self._lo_to_splitting = widgets.Dropdown(
            description="lo_to_splitting:",
            options=["automatic", True, False],
            value="automatic",
        )

        # "Integration method for DOS"
        self._dos_method = widgets.Dropdown(
            description="dos_method:",
            options=["tetra", "gaussian"],
            value="tetra",
        )

        self._stacked_pjdos = widgets.Checkbox(description="Stacked PJDOS", value=True)

        # "Temperature range in K."
        self._temp_range = widgets.FloatRangeSlider(
            description="Temp:",
            value=[0.0, 300.0],
            min=0.0,
            max=1000.0,
            step=0.1,
            readout_format=".1f",
        )

        self._units = widgets.Dropdown(
            description="Energy units:",
            options=["eV", "meV", "Ha", "cm-1", "Thz"],
            value="eV",
        )

        self.children = [
            self._header,
            self._caption,
            self._nqsmall,
            self._ndivsm,
            self._asr,
            self._chneut,
            self._dipdip,
            self._lo_to_splitting,
            self._dos_method,
            self._stacked_pjdos,
            self._temp_range,
            self._units,
        ]
        # traitlets.link((self.model, 'password'), (self._password_text, 'value'))

    @property
    def nqsmall(self):
        return self._nqsmall.value

    @property
    def ndivsm(self):
        return self._ndivsm.value

    @property
    def asr(self):
        return self._asr.value

    @property
    def chneut(self):
        return self._chneut.value

    @property
    def dipdip(self):
        return self._dipdip.value

    @property
    def dos_method(self):
        return self._dos_method.value

    @property
    def lo_to_splitting(self):
        return self._lo_to_splitting.value

    @property
    def stacked_pjdos(self):
        return self._stacked_pjdos.value

    @property
    def temp_range(self):
        return self._temp_range.value

    @property
    def units(self):
        return self._units.value


class PHBandsWidget(widgets.Output):
    def update(self, phbands, phdos, units):
        self.clear_output()
        with self:
            fig = phbands.plotly_with_phdos(phdos, units=units, show=False)
            fig.update_layout(
                title="Phonon band structure and DOS",
            )
            fig.show()


class BrillouinWidget(widgets.HBox):
    def __init__(self):
        super().__init__()

        self._plot = widgets.Output(
            layout={
                "width": "70%",
            }
        )
        self._table = widgets.Output(
            layout={
                "width": "30%",
            }
        )

        self.children = [self._plot, self._table]

    def update(self, phbands):
        self._table.clear_output()
        with self._table:
            display(phbands.qpoints.get_highsym_datataframe())

        self._plot.clear_output()
        with self._plot:
            fig = phbands.qpoints.plotly()
            fig.update_layout(title="Brillouin zone and q-path:")
            fig.show()


class TypeProjectedWidget(widgets.Output):
    def update(self, phdos_file, units, stacked_pjdos):
        self.clear_output()
        with self:
            fig = phdos_file.plotly_pjdos_type(
                units=units, stacked=stacked_pjdos, show=False
            )
            fig.update_layout(title="Type-projected phonon DOS")
            fig.show()


class ThermodynamicWidget(widgets.Output):
    def update(self, phdos, tstart, tstop):
        self.clear_output()
        with self:
            fig = phdos.plotly_harmonic_thermo(
                tstart=tstart, tstop=tstop, num=50, show=False
            )
            fig.update_layout(
                title="Thermodynamic properties in the harmonic approximation"
            )
            fig.show()


class InputWidget(widgets.Output):
    def update(self, inp):
        self.clear_output()
        with self:
            display(widgets.HTML(value="<h2>Anaddb input file:<h2>"))
            display(inp)


class PHBandsPanel(widgets.VBox):
    def __init__(self, ddb):
        super().__init__()

        # filename can be used here
        # ddb = DdbFile.from_file(filename)
        self.ddb = ddb

        self._options = PHBandsOptions()
        self._button = widgets.Button(
            description="Plot Bands and DOS",
            button_style="success",  # 'success', 'info', 'warning', 'danger' or ''
            tooltip="Compute phonon bands and DOSes from DDB and plot the results.",
        )
        self._button.on_click(self.on_click)

        self._phbands = PHBandsWidget()
        self._brillouin = BrillouinWidget()
        self._typeprojected = TypeProjectedWidget()
        self._thermodynamic = ThermodynamicWidget()
        self._input = InputWidget()

        self.children = [
            widgets.VBox([self._options, self._button]),
            widgets.VBox(
                [
                    self._phbands,
                    self._brillouin,
                    self._typeprojected,
                    self._thermodynamic,
                    self._input,
                ]
            ),
        ]

    def update(self):
        """Compute phonon bands and DOSes from DDB and plot the results."""
        options = self._options

        with self.ddb.anaget_phbst_and_phdos_files(
            nqsmall=options.nqsmall,
            ndivsm=options.ndivsm,
            asr=options.asr,
            chneut=options.chneut,
            dipdip=options.dipdip,
            dos_method=options.dos_method,
            lo_to_splitting=options.lo_to_splitting,
            verbose=None,
            mpi_procs=1,
            return_input=True,
        ) as results:

            phbst_file, phdos_file = results
            phbands, phdos = phbst_file.phbands, phdos_file.phdos

            self._phbands.update(phbands, phdos, options.units)
            self._brillouin.update(phbands)
            self._typeprojected.update(phdos_file, options.units, options.stacked_pjdos)
            tstart, tstop = options.temp_range
            self._thermodynamic.update(phdos, tstart, tstop)
            self._input.update(results.input)

    def on_click(self, *args):

        self._button.disabled = True
        self.update()
        self._button.disabled = False
