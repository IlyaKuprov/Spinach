#pragma once

namespace SpinachGUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/** \class  BalloonToolTip
	*  BalloonToolTip class, for tooltips that stay on the top and have better behavior
	* than the normal tooltips as they appear as windows forms
	*/
	public ref class BalloonToolTip : public System::Windows::Forms::Form
	{
	public:
		BalloonToolTip(String ^ ToolTipText)
		{
			InitializeComponent();

			//Text if the ballon
			Ballonlabel->Text=ToolTipText;
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~BalloonToolTip()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  Ballonlabel;
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(BalloonToolTip::typeid));
			this->Ballonlabel = (gcnew System::Windows::Forms::Label());
			this->SuspendLayout();
			// 
			// Ballonlabel
			// 
			this->Ballonlabel->AutoSize = true;
			this->Ballonlabel->BackColor = System::Drawing::Color::Transparent;
			this->Ballonlabel->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->Ballonlabel->Location = System::Drawing::Point(3, 11);
			this->Ballonlabel->Name = L"Ballonlabel";
			this->Ballonlabel->Size = System::Drawing::Size(174, 15);
			this->Ballonlabel->TabIndex = 0;
			this->Ballonlabel->Text = L"Add atom by selecting element";
			// 
			// BalloonToolTip
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
			this->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Stretch;
			this->ClientSize = System::Drawing::Size(208, 35);
			this->Controls->Add(this->Ballonlabel);
			this->FormBorderStyle = System::Windows::Forms::FormBorderStyle::None;
			this->Name = L"BalloonToolTip";
			this->ShowIcon = false;
			this->ShowInTaskbar = false;
			this->StartPosition = System::Windows::Forms::FormStartPosition::Manual;
			this->Text = L"BalloonToolTip";
			this->TransparencyKey = System::Drawing::Color::White;
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	};
}
