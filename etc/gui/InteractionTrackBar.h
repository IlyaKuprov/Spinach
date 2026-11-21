#pragma once
#include <windows.h>
#include <commctrl.h>

namespace SpinachGUI
{
    // Native horizontal slider, with direct selection instead of paging on a track click.
    ref class InteractionTrackBar : public System::Windows::Forms::TrackBar
    {
    protected:
        virtual void WndProc(System::Windows::Forms::Message% message) override
        {
            if (message.Msg == WM_LBUTTONDOWN)
            {
                HWND window = static_cast<HWND>(Handle.ToPointer());
                RECT channel = {}, thumb = {};
                SendMessage(window, TBM_GETCHANNELRECT, 0, reinterpret_cast<LPARAM>(&channel));
                SendMessage(window, TBM_GETTHUMBRECT, 0, reinterpret_cast<LPARAM>(&thumb));
                int x = static_cast<short>(LOWORD(message.LParam.ToInt64()));
                if (x < thumb.left || x >= thumb.right)
                {
                    int halfThumb = (thumb.right - thumb.left) / 2;
                    int start = channel.left + halfThumb;
                    int end = channel.right - halfThumb - 1;
                    if (end > start)
                    {
                        double fraction = System::Math::Min(1.0, System::Math::Max(0.0, (x-start)/double(end-start)));
                        Value = Minimum + static_cast<int>(System::Math::Round(fraction*(Maximum-Minimum)));
                    }
                }
            }
            TrackBar::WndProc(message);
        }
    };
}
